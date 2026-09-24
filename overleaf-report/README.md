# Report source layout

Compile `main.tex`. It contains the document settings and defines the order of the report.

- `chapters/00-abstract.tex`: abstract
- `chapters/01-introduction.tex` through `chapters/12-acknowledgements.tex`: report chapters
- `appendices/`: supporting appendices
- `figures/`: images included by the chapter files
- `references.bib`: bibliography database

To add a chapter, create a numbered `.tex` file in `chapters/` and add one `\input{...}` line in `main.tex`. Labels and references can be used across chapter files without any change.
