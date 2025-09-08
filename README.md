# Pixi Unleashed

This is a replica of the `Elixir-SCO-spatial-omics` repository, but with the addition of a `Pixi` environment for all practicals except practical 5, which is R-based.
The `Pixi` environment is designed to facilitate the execution of Python-based practicals, without the need for switching `Docker` images, which also decreases redundancy of repeating packages = Lightweight & reproducible

The environment also natively supports macOS ARM architecture, making it suitable for users with Apple Silicon devices.

## Usage

To use the `Pixi` environment, follow these steps:
1. Make sure Pixi is installed on your system. You can find installation instructions [here](https://pixi.bio/docs/installation).
2. Clone this repository to your local machine.
3. Navigate to the main directory of the repo .
4. Run the `Jupyter Lab` command to start the Jupyter Lab environment with the `Pixi` environment activated.
5. Open the notebook for the practical you want to work on and that's it!
6. For practical 5, which is R-based, ensure you have the appropriate R environment set up as per the original repository instructions (Docker)

```bash
git clone https://github.com/addityea/ELIXIR-SCO-spatial-omics.git
cd ELIXIR-SCO-spatial-omics
pixi run jupyter-lab
```
