Building the docs locally:

(From this directory.)

1. Run Doxygen
   ```
   doxygen
   ```

2. Run Doxysphinx
   ```
   doxysphinx -v INFO build . ./_build Doxyfile
   ```

3. Run Sphinx
   ```
   sphinx-build -b html . ./_build
   ```

Or, in one line:

```bash
doxygen && doxysphinx -v INFO build . ./_build Doxyfile && sphinx-build -b html . ./_build
```
