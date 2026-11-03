format:
	uv run ruff format

lint:
	uv run ruff check --fix

test:
	uv run pytest

validate: format lint test
