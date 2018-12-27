# splendid

Branch|[![Travis CI logo](pics/TravisCI.png)](https://travis-ci.org)|[![AppVeyor logo](pics/AppVeyor.png)](https://www.appveyor.com)|[![Codecov logo](pics/Codecov.png)](https://www.codecov.io)
---|---|---|---
`master`|[![Build Status](https://travis-ci.org/richelbilderbeek/splendid.svg?branch=master)](https://travis-ci.org/richelbilderbeek/splendid) |[![Build status](https://ci.appveyor.com/api/projects/status/o6htu70cv6ttqw5r/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/splendid/branch/master)| [![codecov.io](https://codecov.io/github/richelbilderbeek/splendid/coverage.svg?branch=master)](https://codecov.io/github/richelbilderbeek/splendid?branch=master)
`develop`|[![Build Status](https://travis-ci.org/richelbilderbeek/splendid.svg?branch=develop)](https://travis-ci.org/richelbilderbeek/splendid) |[![Build status](https://ci.appveyor.com/api/projects/status/o6htu70cv6ttqw5r/branch/develop?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/splendid/branch/develop)| [![codecov.io](https://codecov.io/github/richelbilderbeek/splendid/coverage.svg?branch=develop)](https://codecov.io/github/richelbilderbeek/splendid?branch=develop)
`giovanni`|[![Build Status](https://travis-ci.org/richelbilderbeek/splendid.svg?branch=giovanni)](https://travis-ci.org/richelbilderbeek/splendid) |[![Build status](https://ci.appveyor.com/api/projects/status/o6htu70cv6ttqw5r/branch/giovanni?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/splendid/branch/giovanni)| [![codecov.io](https://codecov.io/github/richelbilderbeek/splendid/coverage.svg?branch=giovanni)](https://codecov.io/github/richelbilderbeek/splendid?branch=giovanni)
`richel`|[![Build Status](https://travis-ci.org/richelbilderbeek/splendid.svg?branch=richel)](https://travis-ci.org/richelbilderbeek/splendid) |[![Build status](https://ci.appveyor.com/api/projects/status/o6htu70cv6ttqw5r/branch/richel?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/splendid/branch/richel)| [![codecov.io](https://codecov.io/github/richelbilderbeek/splendid/coverage.svg?branch=richel)](https://codecov.io/github/richelbilderbeek/splendid?branch=richel)
`pedro`|[![Build Status](https://travis-ci.org/richelbilderbeek/splendid.svg?branch=pedro)](https://travis-ci.org/richelbilderbeek/splendid) |[![Build status](https://ci.appveyor.com/api/projects/status/o6htu70cv6ttqw5r/branch/pedro?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/splendid/branch/pedro)| [![codecov.io](https://codecov.io/github/richelbilderbeek/splendid/coverage.svg?branch=pedro)](https://codecov.io/github/richelbilderbeek/splendid?branch=pedro)

![Didapper](pics/didapper2.jpg)

`splendid` is a SPeciation Likelihood Engine with a DIDapper logo.
(N.B.: `splendid` is a provisional name)

The project is run by G. Laudanno, R.J.C. Bilderbeek and P.M. Santos Neves.

### Goal
Our goal is to make easier to build R packages for likelihood models in macroevolution.
The package is built in a fully modular fashion. In this way the user can build a likelihood package providing:
- one (or more) loglik function(s);
- one (or more) conditioning function(s);
- one (or more) simulation condition(s);
- a function for any event that can occur in simulations;

Once these functions are provided the user should be able to maximize the likelihood and infer the best parameters for any model.
