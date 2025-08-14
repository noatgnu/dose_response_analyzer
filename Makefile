.PHONY: lint format test check-all install dev-install

# Install dependencies
install:
	poetry install

# Install with dev dependencies
dev-install:
	poetry install --with dev

# Format code
format:
	poetry run black .
	poetry run isort .

# Lint code
lint:
	poetry run flake8 .
	poetry run mypy .

# Run tests
test:
	poetry run pytest --cov=dra --cov-report=html --cov-report=term

# Run all checks (format, lint, test)
check-all: format lint test
	@echo "✅ All checks passed!"

# Clean up
clean:
	rm -rf .coverage htmlcov/ .pytest_cache/ .mypy_cache/
	find . -type d -name __pycache__ -delete
	find . -type f -name "*.pyc" -delete
