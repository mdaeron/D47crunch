import tomllib
import pathlib
import datetime

srcpath = next(pathlib.Path("src").glob("*/__init__.py")).parent

with open('pyproject.toml', 'rb') as fid:
	toml = tomllib.load(fid)

export = []
with open(srcpath / '_metadata.py', 'w') as fid:

	align = 16
	for field, content in (
		('version',     toml['project']['version']),
		('description', toml['project']['description']),
		(
			('author', toml['project']['authors'][-1]['name'])
			if len(toml['project']['authors']) == 1 else
			('authors', ' & '.join([', '.join(_['name'] for _ in toml['project']['authors'][:-1]), toml['project']['authors'][-1]['name']]))
		),
		('contact',     toml['project']['authors'][0]['email']),
		('license',     toml['project']['license']),
		('copyright',   toml['project']['authors'][0]['name']),
		('date',        datetime.date.today().isoformat()),
	):
		if content:
			name = f'__{field}__'
			fid.write(f'{name:<{align}s}= "{content}"\n')
			export.append(f'"{name}"')
	fid.write(f'{"__all__":<{align}s}= [{', '.join(export)}]\n')

with open(srcpath / '_metadata.py', 'r') as fid:
	print(fid.read())
