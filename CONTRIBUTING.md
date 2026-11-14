# Contributing to ChargeFW

Bug reports, questions, and pull requests are welcome on
[GitHub](https://github.com/sb-ncbr/chargefw).

## Reporting a problem

Open an issue with the ChargeFW version (`chargefw --version`, or `chargefw.__version__` in Python),
how you installed it, and the command or code you ran. Attach the input file if you can share it. For charges you think are
wrong, name the method and the parameter set.

## Changing ChargeFW

[DEVELOPMENT.md](DEVELOPMENT.md) explains how to build and test ChargeFW, and [AGENTS.md](AGENTS.md)
has the repository's rules for architecture, tests, and documentation. Before you open a pull
request:

- build `gcc-debug` and run its tests, as described in
  [Build and test](DEVELOPMENT.md#build-and-test);
- run `pre-commit run --all-files` and Ruff on the Python files you changed, as described in
  [Formatting](DEVELOPMENT.md#formatting) and [Python checks](DEVELOPMENT.md#python-checks);
- when users will notice the change, update the documentation under [`docs/`](docs/) and add an entry
  under `## Unreleased` in [CHANGELOG.md](CHANGELOG.md).

CI runs the same checks and the wheel tests on every pull request. A pull request can be merged into
`main` only when the formatting checks, Ruff, and the `gcc-debug` tests pass.
