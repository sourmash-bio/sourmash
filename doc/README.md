# Documentation on the docs

We use
[MyST](https://myst-parser.readthedocs.io/en/latest/sphinx/intro.html)
to generate Sphinx doc output from Markdown input.

## Useful tips and tricks:

### Linking internally between sections in the docs

For linking within the sourmash docs, you should use the
[auto-generated header anchors](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#auto-generated-header-anchors)
provided by MyST.

You can generate a list of these for a given document with:
```
myst-anchors -l 3 command-line.md
```
