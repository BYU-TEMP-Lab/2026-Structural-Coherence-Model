from tc_batch_cli import run_many
import batch_plot_configs

if __name__ == "__main__":
    configs = batch_plot_configs.SPEC_CONFIGS
    run_many(configs)
