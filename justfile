default: bayes

bayes:
	cd examples; uv run bayes-example.py

metadata:
	uv run build-metadata.py

doc:
	uv run build_doc.py

diff:
	git diff --stat -- ':!*.png' ':!*.html' ':!*.pdf' ':!*.csv'

publish:
	uv run flit publish
