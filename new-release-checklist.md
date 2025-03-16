# New Release Checklist

1. Bump `__version__`
2. Update `changelog.md`
3. `pixi run docs`
4. `git commit`
5. `git push`
6. Create new release on GitHub
7. `flit publish`