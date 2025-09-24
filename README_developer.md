# Developer guide

> [!TIP]
> Always work from the root directory of the repository.

## Setup

Always use the latest version of Python when developing this project.

This project uses [uv](https://docs.astral.sh/uv/) for dependency and project management.
To install uv, follow the [install instructions](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer).

A few of the most useful uv commands are shown below.
See the [uv project guide](https://docs.astral.sh/uv/guides/projects) for more information about working on projects with uv.

- `uv add` -- Add dependencies to the project
- `uv remove` -- Remove dependencies from the project
- `uv sync` -- Update the project's environment
- `uv run` -- Run a command or script

> [!TIP]
> Instead of `pip install <package>`, use `uv add <package>`.

1.  **Clone** this repo and run `uv sync --dev` to install the project and development dependencies.

> [!NOTE]
> Note that uv creates a virtual environment located at `./venv` to store the project dependencies.
> Explicit activation of this virtual environment is _not necessary_ thanks to the `uv run` command.

## Scripts

To run your code, or a command that is aware of your project environment, use the following command:

- `uv run` -- Run a command or script

See the [script documentation](https://docs.astral.sh/uv/guides/scripts/) for more information.

### Code quality

Use [Ruff](https://docs.astral.sh/ruff/), [mypy](https://www.mypy-lang.org/), and [pytest](https://docs.pytest.org/en/stable/)to respectively lint/format, typecheck, and test your code.
These tools are already specified as development dependencies in this template's `pyproject.toml`, and are automatically installed when you run `uv sync --dev`.

> [!NOTE]
> Scripts to check/enforce code quality are supplied in `./scripts/`; run them with `uv run scripts/{lint,format,test}`.
> You can also run these tools or others individually as needed.

There are many options available to customize and configure these tools.
See `pyproject.toml` for more information.

## Branching Strategy

- The `master` branch always has deployable code and never receives direct commits.
- Releases are created as tags on the `master` branch.
- All changes are made in small chunks on specific feature branches and incorporated to `master` via squashed pull request.

## Publishing a new Version

Please see our Group's Wiki if you are looking to publish a new version.

## Updating Documentation

Documentation is automatically created by our flask site. Comment blocks above functions are parsed into html pages.
The format of the comment blocks is provided in the **init**.py script in the src folder. TIP: If you want to exclude a function's comment block
from being displayed on the flask site, use single quotes ''' instead of double quotes """.

## Github Actions

We have two github actions setup. One is used to publish the code to Pypi and the other is used to run
our test scripts across multiple versions. These actions are defined in the .github/workflow folder. The associated tests
that are run with each actions are in the Tests folder titled test*github*\*.py. SEE ALSO : README.md in Tests.
