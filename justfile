default: start

start:
	open src/D47crunch/__init__.py

bayes:
	cd examples; uv run bayes-demo.py

metadata:
	uv run build-metadata.py

doc:
	uv run build_doc.py

diff:
	git diff --stat -- ':!*.png' ':!*.html' ':!*.pdf' ':!*.csv'

publish:
	uv build && uv publish

test:
	uv run pytest tests -s
