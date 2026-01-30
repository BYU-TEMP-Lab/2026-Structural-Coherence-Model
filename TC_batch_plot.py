from TC_batch_cli import run_many
import TC_batch_plot_configs

if __name__ == "__main__":
    configs = TC_batch_plot_configs.SPEC_CONFIGS
    run_many(configs)
