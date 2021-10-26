% Rd macro for simplifying documentation writing

\newcommand{\fun}{\code{\link[=#1]{#1()}}}
%E.g., \fun{mean} is the same as \code{link[=mean]{mean}}

\newcommand{\pkgfun}{\code{\link[#1:#2]{#2()}}}

%E.g., \pkgfun{sandwich}{vcovCL} is the same as \code{\link[sandwich:vcovCL]{vcovCL}}

\newcommand{\pkgfun2}{\code{\link[#1:#2]{#3()}}}

%E.g., \pkgfun2{sandwich}{vcovCL}{meatCL} is the same as \code{\link[sandwich:vcovCL]{meatCL}}