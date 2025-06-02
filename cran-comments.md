## Comments on first submission notes

These are comments on the notes given by a CRAN member on the first tried submission of ```deseats``` (November 03, 2023).

* Please always explain all acronyms in the description text. e.g.: TS.

The acronym TS in the package title is now explained in the description text.

* It seems as you have too many blank spaces in your description text. Probably because line breaks count as blank spaces too. Please omit this.

Unnecessary blank spaces at the end of lines in the description text have been removed, so that the line breaks work as individual blank spaces and not double blank spaces.

* Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
  Missing Rd-tags:
     animate.Rd: \value
     expo.Rd: \value
     gain.Rd: \value
     order_poly.Rd: \value
     trend.Rd: \value
     
These Rd-\value-tags were added.

* \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
  Does not seem necessary.
  Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}.

\dontrun{} was either removed throughout or replaced with donttest{}.

* You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object.
  Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console.
  (except for print, summary, interactive functions).

```cat()``` in ```bwidth_confint()``` was removed. All other calls to ```cat()``` in the ```deseats``` package are either in ```show()``` methods or in an interactive context.

* Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies.
  Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir().
  
The example in ```read_ts()``` was fixed accordingly, so that it now writes to
```tempdir()```.

* You are setting options(warn=-1) in your function. This is not allowed.  e.g.: inst/shiny-apps/lwr/app.R

The line ```options(warn=-1)``` was removed from the app file.
