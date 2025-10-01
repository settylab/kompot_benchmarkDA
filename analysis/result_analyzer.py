#!/usr/bin/env python3
"""
BenchmarkDA Result Analyzer
============================

Modern analysis framework for BenchmarkDA results that works with the new infrastructure.
Supports flexible data discovery, comprehensive metrics, and modular plotting functions.

Key Features:
- Automatic result discovery from new hierarchical structure
- Configurable dataset/method/embedding filtering
- Comprehensive metrics: AUROC, AUPRC, correlation, ranking metrics
- Modular plotting functions for different visualizations
- Ground truth integration and validation
"""

import os
import glob
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.metrics import (
    auc, roc_auc_score, precision_recall_curve,
    average_precision_score, mean_absolute_error,
    mean_squared_error, r2_score
)
from scipy.stats import pearsonr, spearmanr, kendalltau, rankdata

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from tqdm.auto import tqdm

# Suppress warnings
warnings.filterwarnings('ignore')

@dataclass
class ResultConfig:
    """Configuration for result analysis."""
    # Paths
    benchmark_root: str = "benchmark"
    data_root: str = "data"

    # Filters
    datasets: Optional[List[str]] = None
    methods: Optional[List[str]] = None
    embeddings: Optional[List[str]] = None
    populations: Optional[List[str]] = None
    enrichments: Optional[List[float]] = None
    seeds: Optional[List[int]] = None
    batch_effects: Optional[List[float]] = None

    # Metrics to compute
    metrics: List[str] = field(default_factory=lambda: [
        'auroc', 'auroc_da', 'auroc_neg', 'auroc_pos',
        'auprc', 'auprc_da', 'auprc_neg', 'auprc_pos',
        'mae', 'mse', 'r2', 'pearson', 'spearman', 'kendall'
    ])

    # Analysis settings
    balance_direction: str = "No"  # Balance direction for job IDs
    iteration: int = 0  # Iteration number

    def __post_init__(self):
        """Set default values based on configuration."""
        if self.datasets is None:
            # Get datasets from config
            try:
                import sys
                from pathlib import Path
                sys.path.insert(0, str(Path.cwd()))
                from config.dataset_config import DATASET_CONFIGS, SEEDS, ENRICHMENT_VALUES
                self.datasets = list(DATASET_CONFIGS.keys())
                if self.seeds is None:
                    self.seeds = SEEDS
                if self.enrichments is None:
                    self.enrichments = ENRICHMENT_VALUES
            except ImportError:
                # Fallback to hardcoded defaults
                self.datasets = ['linear', 'branch', 'cluster', 'covid19-pbmc', 'bcr-xl', 'levine32', 'pancreas']
                if self.seeds is None:
                    self.seeds = [43, 44, 45]
                if self.enrichments is None:
                    self.enrichments = [0.75, 0.85, 0.95]

        if self.embeddings is None:
            self.embeddings = ['dm', 'pca']
        if self.batch_effects is None:
            self.batch_effects = [0.0, 0.75, 1.0, 1.25, 1.5]

class ResultAnalyzer:
    """Main analyzer class for BenchmarkDA results."""

    def __init__(self, config: ResultConfig):
        self.config = config
        self.project_root = Path.cwd()
        self._load_dataset_info()

    def _load_dataset_info(self):
        """Load dataset configuration information."""
        try:
            import sys
            sys.path.append(str(self.project_root))
            from config.dataset_config import DATASET_CONFIGS
            self.dataset_configs = DATASET_CONFIGS
        except ImportError:
            print("Warning: Could not load dataset configuration. Using default population lists.")
            self.dataset_configs = self._get_default_populations()

    def _get_default_populations(self) -> Dict:
        """Fallback population definitions."""
        return {
            'linear': {'pops': [f'M{i}' for i in range(1, 8)]},
            'branch': {'pops': [f'M{i}' for i in range(1, 9)]},
            'cluster': {'pops': [f'M{i}' for i in range(1, 4)]},
            'covid19-pbmc': {'pops': ['B', 'RBC', 'PB', 'CD14_Monocyte', 'CD8_T', 'CD4_T', 'Platelet', 'NK', 'Granulocyte', 'CD16_Monocyte', 'gd_T', 'pDC', 'DC']},
            'bcr-xl': {'pops': ['CD4_T-cells', 'NK_cells', 'CD8_T-cells', 'monocytes', 'B-cells_IgM-', 'DC', 'B-cells_IgM+', 'surface-']},
            'levine32': {'pops': ['pDCs', 'CD4_T_cells', 'CD8_T_cells', 'Pre_B_cells', 'Mature_B_cells', 'Monocytes', 'Basophils']},
            'pancreas': {'pops': ['delta_cell', 'alpha_cell', 'gamma_cell', 'acinar_cell', 'beta_cell', 'ductal_cell', 'epsilon_cell']}
        }

    def discover_results(self) -> Dict:
        """Discover available results in the new infrastructure."""
        results = {}
        benchmark_root = self.project_root / self.config.benchmark_root

        if not benchmark_root.exists():
            print(f"Benchmark directory not found: {benchmark_root}")
            return results

        print("🔍 Discovering results...")

        for dataset in tqdm(self.config.datasets, desc="Datasets"):
            # Check both synthetic and real data paths
            for data_type in ['synthetic', 'real']:
                dataset_path = benchmark_root / data_type / dataset
                if not dataset_path.exists():
                    continue

                # Find all job directories
                job_dirs = [d for d in dataset_path.iterdir() if d.is_dir()]

                for job_dir in job_dirs:
                    # Parse job directory name: dataset-pop-enr-seed-batch-balance-embedding
                    job_parts = job_dir.name.split('-')
                    if len(job_parts) < 7:
                        continue

                    job_dataset, pop, enr, seed, batch, balance, embedding = job_parts[:7]

                    # Apply filters
                    if self.config.populations and pop not in self.config.populations:
                        continue
                    if self.config.embeddings and embedding not in self.config.embeddings:
                        continue
                    if self.config.enrichments and float(enr) not in self.config.enrichments:
                        continue
                    if self.config.seeds and int(seed) not in self.config.seeds:
                        continue
                    if self.config.batch_effects and float(batch) not in self.config.batch_effects:
                        continue

                    # Look for iteration directory
                    iter_dir = job_dir / f"iteration_{self.config.iteration}"
                    if not iter_dir.exists():
                        continue

                    # Find result files
                    result_files = list(iter_dir.glob("*.csv"))
                    if not result_files:
                        continue

                    # Store results
                    job_key = f"{dataset}_{pop}_{enr}_{seed}_{batch}_{balance}_{embedding}"

                    if dataset not in results:
                        results[dataset] = {}
                    if job_key not in results[dataset]:
                        results[dataset][job_key] = {
                            'data_type': data_type,
                            'population': pop,
                            'enrichment': float(enr),
                            'seed': int(seed),
                            'batch_effect': float(batch),
                            'embedding': embedding,
                            'job_dir': job_dir,
                            'iter_dir': iter_dir,
                            'result_files': result_files,
                            'methods': {}
                        }

                    # Parse result files by method
                    for result_file in result_files:
                        if 'DAresults' in result_file.name:
                            method = result_file.name.split('.')[-2]  # Extract method name
                            if self.config.methods is None or method in self.config.methods:
                                results[dataset][job_key]['methods'][method] = result_file

        total_jobs = sum(len(jobs) for jobs in results.values())
        print(f"📊 Found {total_jobs} result combinations across {len(results)} datasets")
        return results

    def load_ground_truth(self, dataset: str, population: str, enrichment: float,
                         seed: int, batch_effect: float = 0.0) -> Optional[pd.DataFrame]:
        """Load ground truth data for a specific parameter combination."""

        # Try to find ground truth in the data generation output
        data_type = 'synthetic' if dataset in ['linear', 'branch', 'cluster'] else 'real'

        # Construct unified job ID for label data
        unified_jobid = f"{dataset}-{population}-{enrichment}-{seed}-{batch_effect}-{self.config.balance_direction}"
        data_dir = self.project_root / self.config.data_root / data_type / dataset / unified_jobid

        # Look for coldata file
        coldata_files = list(data_dir.glob("*.coldata.csv")) if data_dir.exists() else []

        if coldata_files:
            coldata_file = coldata_files[0]
            try:
                coldata = pd.read_csv(coldata_file, index_col=0)

                # Extract ground truth information
                if all(col in coldata.columns for col in ['true_labels', 'Condition1_prob', 'Condition2_prob']):
                    ground_truth = pd.DataFrame(index=coldata.index)
                    ground_truth['true_labels'] = coldata['true_labels']
                    ground_truth['condition1_prob'] = coldata['Condition1_prob']
                    ground_truth['condition2_prob'] = coldata['Condition2_prob']
                    ground_truth['synth_labels'] = coldata.get('synth_labels', 'Unknown')

                    # Calculate log fold change
                    ground_truth['lfc_truth'] = np.log(ground_truth['condition2_prob']) - np.log(ground_truth['condition1_prob'])

                    # Create DA labels (any change)
                    ground_truth['da_labels'] = ground_truth['true_labels'].apply(
                        lambda x: 'DA' if x in ['NegLFC', 'PosLFC'] else 'NotDA'
                    )

                    return ground_truth
            except Exception as e:
                print(f"Warning: Could not load ground truth for {unified_jobid}: {e}")

        return None

    def load_method_results(self, result_file: Path) -> Optional[pd.DataFrame]:
        """Load results from a specific method result file."""
        try:
            return pd.read_csv(result_file, index_col=0)
        except Exception as e:
            print(f"Warning: Could not load results from {result_file}: {e}")
            return None

    def compute_metrics(self, ground_truth: pd.DataFrame,
                       predictions: pd.DataFrame, method: str) -> Dict[str, float]:
        """Compute all specified metrics for a method's predictions."""
        metrics = {}

        if predictions.empty or method not in predictions.columns:
            return {metric: np.nan for metric in self.config.metrics}

        y_pred = predictions[method].values
        y_true_lfc = ground_truth['lfc_truth'].values

        # Get label arrays
        true_labels = ground_truth['true_labels'].values
        da_labels = ground_truth['da_labels'].values

        true_neg = (true_labels == 'NegLFC')
        true_pos = (true_labels == 'PosLFC')
        true_da = (da_labels == 'DA')

        # Classification metrics
        if 'auroc' in self.config.metrics:
            auroc_neg = self._safe_auroc(true_neg, -y_pred)
            auroc_pos = self._safe_auroc(true_pos, y_pred)
            metrics['auroc_neg'] = auroc_neg
            metrics['auroc_pos'] = auroc_pos
            metrics['auroc'] = np.nanmean([auroc_neg, auroc_pos])

        if 'auroc_da' in self.config.metrics:
            metrics['auroc_da'] = self._safe_auroc(true_da, np.abs(y_pred))

        if 'auprc' in self.config.metrics:
            auprc_neg = self._safe_auprc(true_neg, -y_pred)
            auprc_pos = self._safe_auprc(true_pos, y_pred)
            metrics['auprc_neg'] = auprc_neg
            metrics['auprc_pos'] = auprc_pos
            metrics['auprc'] = np.nanmean([auprc_neg, auprc_pos])

        if 'auprc_da' in self.config.metrics:
            metrics['auprc_da'] = self._safe_auprc(true_da, np.abs(y_pred))

        # Regression metrics
        if 'mae' in self.config.metrics:
            metrics['mae'] = mean_absolute_error(y_true_lfc, y_pred)

        if 'mse' in self.config.metrics:
            metrics['mse'] = mean_squared_error(y_true_lfc, y_pred)

        if 'r2' in self.config.metrics:
            r2 = r2_score(y_true_lfc, y_pred)
            metrics['r2'] = max(0, r2)  # Clip negative R²

        # Correlation metrics
        if 'pearson' in self.config.metrics:
            corr, _ = pearsonr(y_true_lfc, y_pred)
            metrics['pearson'] = corr

        if 'spearman' in self.config.metrics:
            corr, _ = spearmanr(y_true_lfc, y_pred)
            metrics['spearman'] = corr

        if 'kendall' in self.config.metrics:
            corr, _ = kendalltau(y_true_lfc, y_pred)
            metrics['kendall'] = corr

        return metrics

    def _safe_auroc(self, y_true: np.ndarray, y_scores: np.ndarray) -> float:
        """Safely compute AUROC, handling edge cases."""
        try:
            if np.sum(y_true) == 0 or np.sum(y_true) == len(y_true):
                return np.nan  # No positive or no negative samples
            return roc_auc_score(y_true, y_scores)
        except:
            return np.nan

    def _safe_auprc(self, y_true: np.ndarray, y_scores: np.ndarray) -> float:
        """Safely compute AUPRC, handling edge cases."""
        try:
            if np.sum(y_true) == 0:
                return np.nan  # No positive samples
            return average_precision_score(y_true, y_scores)
        except:
            return np.nan

    def analyze_all_results(self) -> pd.DataFrame:
        """Analyze all discovered results and return comprehensive metrics DataFrame."""
        results = self.discover_results()

        if not results:
            print("No results found!")
            return pd.DataFrame()

        all_metrics = []

        print("🧮 Computing metrics...")

        for dataset in tqdm(results.keys(), desc="Datasets"):
            dataset_results = results[dataset]

            for job_key, job_info in tqdm(dataset_results.items(), desc=f"{dataset} jobs", leave=False):
                # Load ground truth
                ground_truth = self.load_ground_truth(
                    dataset, job_info['population'], job_info['enrichment'],
                    job_info['seed'], job_info['batch_effect']
                )

                if ground_truth is None:
                    continue

                # Process each method
                for method, result_file in job_info['methods'].items():
                    method_results = self.load_method_results(result_file)

                    if method_results is None:
                        continue

                    # Compute metrics
                    metrics = self.compute_metrics(ground_truth, method_results, method)

                    # Create record
                    record = {
                        'dataset': dataset,
                        'data_type': job_info['data_type'],
                        'population': job_info['population'],
                        'enrichment': job_info['enrichment'],
                        'seed': job_info['seed'],
                        'batch_effect': job_info['batch_effect'],
                        'embedding': job_info['embedding'],
                        'method': method,
                        'job_key': job_key,
                        **metrics
                    }
                    all_metrics.append(record)

        metrics_df = pd.DataFrame(all_metrics)
        print(f"📈 Computed metrics for {len(metrics_df)} method-job combinations")

        return metrics_df

    def plot_performance_comparison(self, metrics_df: pd.DataFrame,
                                  metric: str = 'auroc',
                                  group_by: str = 'enrichment',
                                  datasets: Optional[List[str]] = None,
                                  methods: Optional[List[str]] = None,
                                  embeddings: Optional[List[str]] = None,
                                  figsize: Tuple[int, int] = (15, 5)) -> plt.Figure:
        """Create performance comparison plots."""

        # Filter data
        plot_df = metrics_df.copy()
        if datasets:
            plot_df = plot_df[plot_df['dataset'].isin(datasets)]
        if methods:
            plot_df = plot_df[plot_df['method'].isin(methods)]
        if embeddings:
            plot_df = plot_df[plot_df['embedding'].isin(embeddings)]

        if plot_df.empty:
            print("No data to plot after filtering")
            return None

        # Get unique datasets for subplots
        unique_datasets = sorted(plot_df['dataset'].unique())
        n_datasets = len(unique_datasets)

        fig, axes = plt.subplots(1, n_datasets, figsize=(figsize[0] * n_datasets, figsize[1]))
        if n_datasets == 1:
            axes = [axes]

        # Create color palette for methods
        unique_methods = sorted(plot_df['method'].unique())
        colors = sns.color_palette("husl", len(unique_methods))
        method_colors = dict(zip(unique_methods, colors))

        for idx, dataset in enumerate(unique_datasets):
            ax = axes[idx]
            dataset_df = plot_df[plot_df['dataset'] == dataset]

            sns.boxplot(data=dataset_df, x=group_by, y=metric,
                       hue='method', ax=ax, palette=method_colors)

            ax.set_title(f'{dataset}')
            ax.set_xlabel(group_by.capitalize())
            ax.set_ylabel(metric.upper())

            if idx < n_datasets - 1:
                ax.get_legend().remove()

        plt.tight_layout()
        return fig

    def plot_method_ranking(self, metrics_df: pd.DataFrame,
                           metric: str = 'auroc',
                           datasets: Optional[List[str]] = None,
                           figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
        """Plot method ranking across datasets."""

        # Filter data
        plot_df = metrics_df.copy()
        if datasets:
            plot_df = plot_df[plot_df['dataset'].isin(datasets)]

        # Compute mean performance per method per dataset
        summary = plot_df.groupby(['dataset', 'method'])[metric].agg(['mean', 'std']).reset_index()

        # Pivot for heatmap
        heatmap_data = summary.pivot(index='method', columns='dataset', values='mean')

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='RdYlBu_r',
                   center=heatmap_data.median().median(), ax=ax)

        ax.set_title(f'Method Performance ({metric.upper()}) Across Datasets')
        ax.set_xlabel('Dataset')
        ax.set_ylabel('Method')

        plt.tight_layout()
        return fig

    def generate_summary_report(self, metrics_df: pd.DataFrame) -> str:
        """Generate a comprehensive summary report."""

        report = ["🎯 BenchmarkDA Analysis Summary", "=" * 50, ""]

        # Basic statistics
        n_datasets = metrics_df['dataset'].nunique()
        n_methods = metrics_df['method'].nunique()
        n_jobs = len(metrics_df)

        report.extend([
            f"📊 **Overview:**",
            f"   • Datasets analyzed: {n_datasets}",
            f"   • Methods evaluated: {n_methods}",
            f"   • Total comparisons: {n_jobs:,}",
            ""
        ])

        # Dataset breakdown
        dataset_counts = metrics_df['dataset'].value_counts().sort_index()
        report.extend([
            f"📈 **Results by Dataset:**"
        ])
        for dataset, count in dataset_counts.items():
            report.append(f"   • {dataset}: {count:,} method evaluations")
        report.append("")

        # Method breakdown
        method_counts = metrics_df['method'].value_counts().sort_index()
        report.extend([
            f"🔧 **Methods Evaluated:**"
        ])
        for method, count in method_counts.items():
            report.append(f"   • {method}: {count:,} evaluations")
        report.append("")

        # Top performers by metric
        for metric in ['auroc', 'auprc', 'pearson']:
            if metric in metrics_df.columns:
                top_methods = (metrics_df.groupby('method')[metric]
                              .mean()
                              .sort_values(ascending=False)
                              .head(3))

                report.extend([
                    f"🏆 **Top 3 Methods by {metric.upper()}:**"
                ])
                for i, (method, score) in enumerate(top_methods.items(), 1):
                    report.append(f"   {i}. {method}: {score:.3f}")
                report.append("")

        return "\n".join(report)


def create_default_config(**kwargs) -> ResultConfig:
    """Create a default configuration with optional overrides."""
    return ResultConfig(**kwargs)


# Example usage and testing functions
def main():
    """Command-line interface for the ResultAnalyzer."""
    import argparse

    parser = argparse.ArgumentParser(
        description="BenchmarkDA Result Analyzer - Comprehensive analysis of differential abundance benchmarking results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick summary of all results
  python result_analyzer.py summary

  # Focus on specific datasets and methods
  python result_analyzer.py compare --datasets linear,cluster --methods mellon,meld

  # Generate detailed report
  python result_analyzer.py report --metrics auroc,auprc
        """
    )

    parser.add_argument('action', nargs='?', default='summary',
                       choices=['summary', 'compare', 'report'],
                       help='Action to perform (default: summary)')

    parser.add_argument('--datasets', type=str,
                       help='Comma-separated list of datasets to analyze')

    parser.add_argument('--methods', type=str,
                       help='Comma-separated list of methods to analyze')

    parser.add_argument('--metrics', type=str, default='auroc,auprc,pearson,spearman',
                       help='Comma-separated list of metrics to compute')

    parser.add_argument('--output', type=str, default='analysis_results.csv',
                       help='Output file for results (default: analysis_results.csv)')

    args = parser.parse_args()

    # Parse comma-separated lists
    datasets = args.datasets.split(',') if args.datasets else None
    methods = args.methods.split(',') if args.methods else None
    metrics = args.metrics.split(',') if args.metrics else ['auroc', 'auprc', 'pearson', 'spearman']

    # Create configuration
    config = create_default_config(
        datasets=datasets,
        methods=methods,
        metrics=metrics
    )

    # Create analyzer
    analyzer = ResultAnalyzer(config)

    # Analyze all results
    print("Starting comprehensive analysis...")
    metrics_df = analyzer.analyze_all_results()

    if metrics_df.empty:
        print("No results found. Make sure benchmarking has been completed.")
        return

    # Handle different actions
    if args.action == 'summary':
        print("\n" + analyzer.generate_summary_report(metrics_df))

    elif args.action == 'compare':
        print("\n" + analyzer.generate_summary_report(metrics_df))
        print("\n🎯 Method Comparison Results:")
        for metric in metrics:
            if metric in metrics_df.columns:
                print(f"\n📊 {metric.upper()} Performance:")
                method_scores = metrics_df.groupby('method')[metric].mean().sort_values(ascending=False)
                for i, (method, score) in enumerate(method_scores.head(5).items(), 1):
                    print(f"  {i}. {method}: {score:.3f}")

    elif args.action == 'report':
        print("\n" + analyzer.generate_summary_report(metrics_df))
        print("\n📈 Generating plots...")
        try:
            import matplotlib.pyplot as plt
            # Performance comparison
            fig1 = analyzer.plot_performance_comparison(
                metrics_df, metric='auroc', group_by='enrichment'
            )
            if fig1:
                fig1.suptitle('AUROC Performance by Enrichment Level')
                plt.show()

            # Method ranking
            fig2 = analyzer.plot_method_ranking(metrics_df, metric='auroc')
            if fig2:
                plt.show()
        except ImportError:
            print("⚠️ Matplotlib not available, skipping plots")

    # Save results
    output_path = Path(args.output)
    metrics_df.to_csv(output_path, index=False)
    print(f"\n💾 Results saved to: {output_path}")


if __name__ == "__main__":
    main()