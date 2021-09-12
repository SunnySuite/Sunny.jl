# Getting started with Julia for users of  _Sunny_

## Julia basics

The first step is to install [Julia](https://julialang.org/). Download the latest stable release (1.6 as of Sept 2021) appropriate to your system from the [downloads page](https://julialang.org/downloads/). Inside the package is a Julia executable `julia` which launches the Julia terminal. Follow the [platform specific instructions](https://julialang.org/downloads/platform/) to make the path to `julia` known to your system.

Executing the command `julia` should bring you to a prompt like this:

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0 (2021-03-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
```

This is called the Julia REPL, short for read-eval-print-loop. You can run code here, for example
```julia
julia> exp(π * im)
-1.0 - 1.2246467991473532e-16im
```
which verifies Euler's identity `e^{π i} = -1` to order `10^-16`.

There are many resources online to learn Julia. To see some examples, there is [_Learn X in Y minutes, where X=Julia_](https://learnxinyminutes.com/docs/julia/) and the [Julia Cheat-sheet](https://juliadocs.github.io/Julia-Cheat-Sheet/). There is also the [official Julia documentation](https://docs.julialang.org/), but you don't need to be a Julia expert get started with Sunny.

Julia has built-in help. From the REPL, enter `?` to enter the help mode, and then enter a function or type name to get documentation. For example,
```
help?> exp
search: exp exp2 Expr expm1 exp10 export exponent expanduser ExponentialBackOff ldexp frexp nextpow

  exp(x)

  Compute the natural base exponential of x, in other words e^x.
  ...
```

Math symbols can be accessed with Latex-like notation. For example, typing `\pi` and then pressing Tab will produce `π`. If you see a unicode symbol that you don't know how to type, you can copy-paste it into the help prompt.

```
help?> ℋ
"ℋ" can be typed by \scrH<tab>
```

Although Julia and Matlab have similar surface syntax, there are important differences between the two languages. One crucial thing to know is that simple "for loops" over scalars are _encouraged_ in Julia, and can be used to generate code with C++ or Fortran-like speed. Vectorized style is supported in Julia, but not necessary. Like Matlab and Fortran (but unlike Python or C++), Julia uses 1-based indexing by default, so `my_array[1]` is the first element of `my_array`.

## The built-in Julia package manager

You may be familiar with package managers like `pip` and `anaconda` for Python. In Julia this functionality is built in.

At the Julia REPL, press the key `]` to enter the "package" mode. You will see that the prompt text changes:

```
(@v1.6) pkg> 
```

The text `@v1.6` indicates that you are now controlling the package "environment" that is used by default when running Julia version 1.6. Typing `help` will give a list of available commands. 

New packages can be added with the `add` command. For example,
```
pkg> add LinearAlgebra
pkg> add StaticArrays
pkg> add OffsetArrays
```
adds three packages that will be useful for Sunny.

There are a lot of interesting Julia packages available. You can browse the _registered_ packages at [juliapackages.com](https://juliapackages.com/), along with their documentation.

Sunny is not yet a registered package, but you can add it directly from its Github URL,
```
pkg> add https://github.com/MagSims/FastDipole.git
```

[[**Note**: This will start working once we (hopefully soon) make Sunny a public repository.]]


On a Unix-based system such as Mac or Linux, Julia stores package information inside the hidden directory `~/.julia`, where "`~`" represents the user home directory. Browse around to see what files Julia created. The file `~/.julia/environments/v1.6/Project.toml` lists all the installed packages for the default Julia 1.6 environment. Another way to get similar information is to enter

```
pkg> status
```

Packages can be removed with the `remove` command, or updated with the `update` command.

To leave the package manager and return to the Julia REPL, press the Backspace key.


## Plotting

The plotting ecosystem in Julia is still evolving rapidly. For simple tasks, the easiest place to start is [Plots](http://docs.juliaplots.org/latest/). Install it with

```
pkg> add Plots
```

Next, return to the REPL by pressing the Backspace key and then enter
```julia
using Plots
# Build a linear array of numbers between [0.1, 20] with step size 0.1
xs = 0.1:0.1:20
# The 'dot' (broadcast) syntax applies the function to all array elements one-by-one
plot(xs, sin.(xs) ./ xs) 
```
This should eventually bring up a plotting window. The "time to first plot" should be of order 20 seconds.

An exciting plotting package, [Makie](https://github.com/JuliaPlots/Makie.jl), is currently under heavy development. It supports interactive 3D graphics and customizable GUI elements like buttons and sliders. Install with
```
pkg> add GLMakie
pkg> test GLMakie
```
There are still rough edges in Makie, so it's a good idea to check that these tests pass. As of Sept. 2021, one annoyance is that it time-to-first-plot is about a minute. Important progress will be made in [this refactoring](https://github.com/JuliaPlots/Makie.jl/pull/1085). Beware also some [longstanding bugs on Mac](https://github.com/JuliaPlots/Makie.jl/milestone/4).

For making publication quality plots, it's hard to beat Matplotlib from Python. Some people use Makie. Other Julia options include [PGFPlotsX](https://github.com/KristofferC/PGFPlotsX.jl) (which generates Latex directly) and [VegaLite](https://github.com/queryverse/VegaLite.jl) (a Julia interface to the Javascript library [Vega-Lite](https://vega.github.io/vega-lite/)). There also exists a [Julia interface to Matplotlib](https://github.com/JuliaPy/PyPlot.jl), but it is probably easier to just use Matplot lib directly from Python.

## Notebook environments

Besides the REPL, another way to interact with Julia is through the Jupyter notebook interface, which may be familiar to Python users. One can integrate Julia with the Jupyter notebook system through the [IJulia](https://github.com/JuliaLang/IJulia.jl) package. Install it like usual from the package manager,

```
pkg> add IJulia
```

Launching a Julia/Jupyter notebook should be as simple as running

```julia
using IJulia
notebook()
```

If Jupyter has already been installed on your system, then IJulia will likely find and use it. Otherwise, IJulia will try to install Jupyter from scratch.

*Caution*: To install Jupyter from scratch, the IJulia package will download the entire Anaconda distribution, which requires gigabytes of storage. If an Anaconda installation already exists on your system, it is probably good to have IJulia find it.

Once all the above packages are installed, you will be ready to run and modify the Sunny examples.

# More advanced features for code development

This section describes some more advanced setup that will be useful for those ready to dig deeper into Sunny.

## Developing the Sunny source code with Git

Previously, we included Sunny into the Julia environment using the package `add` command. For advanced users, it will be better to use the `develop` command, which allows modifying the Sunny source code. If Sunny has been installed with `add`, remove it using

```
pkg> remove FastDipole
```

Next run
```
pkg> develop https://github.com/MagSims/FastDipole.git
```

On a Unix system, this will download the Sunny repository into the local directory `~/.julia/dev/FastDipole`. You can freely modify any of the source code in the `~/.julia/dev/` directory, and your changes will be picked up by Julia.

The [Git version control system](https://git-scm.com/) makes it possible to fearlessly experiment with code modifications. A full introduction to Git is beyond the scope of this document, but here are some basics. Try modifying any file inside the `FastDipole/` directory. Next, open a terminal in the `FastDipole/` directory and type `git status`. You should see the name of the file that was changed. Type `git diff` to see the specific changes made. You can revert these changes with the command `git checkout <filename>`. Other useful commands include `git add <filenames>` and `git commit`, which will enter changes into the database (repository) of tracked changes (commits). `git log` will show a history of commits. The commands `git pull` and `git push` will download and upload, respectively, from the "origin" repository (in this case, the one hosted on Github). To have new commits pushed into the team Github repo, we will use a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/) workflow. All proposed changes will be submitted as pull requests, where they will be automatically tested and granted code review.

Instead of `git` terminal commands, it's often more convenient to use a graphical user interface. For Julia development, I would highly recommend the VSCode editor, described below. VSCode has some built-in support for Git. It may also be useful to have a dedicated Git GUI such as [SmartGit](https://www.syntevo.com/smartgit/), which free for open source use.

## Revise-based workflow

An extremely important package in the Julia ecosystem is [Revise](https://timholy.github.io/Revise.jl/stable/),
```
pkg> add Revise
```

Revise allows a running Julia process to automatically pick up changes to package source code (i.e., to dynamically replace function behaviors) without requiring a restart. This is very convenient for rapid iteration -- one can avoid most of the wait times associated with reloading (and recompiling) packages.

To enable Revise by default, create (or extend) the file `~/.julia/config/startup.jl`. The content of this file is automatically executed at the start of every Julia REPL sesion. Include the following lines:
```
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
```

To enable Revise support in Jupyter notebooks, similar lines should be added to `~/.julia/config/startup_ijulia.jl`,
```
try
    @eval using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
```

Revise works great for redefining function definitions. A limitation of Revise is that it does _not_ support redefining `struct` datatypes. In such cases, the Julia process must (unfortunately) be restarted.

## The Julia extension in VSCode

[VSCode](https://code.visualstudio.com/) is a powerful text editor, and hosts the _de facto_ official [Julia development environment](https://www.julia-vscode.org/). Note that VSCode is [open source](https://github.com/microsoft/vscode), and is entirely unrelated to the commercial Visual Studio product. (It seems Microsoft was deliberately trying create confusion with this naming scheme.)

A full introduction to VSCode is outside the scope of this document, but here we can provide a few helpful tips. Much of the power of VSCode comes from extensions. You can download and install these easily through the Extensions panel, accessible by clicking the appropriate icon on the left of the VSCode window (alternatively, by selecting the `View -> Extensions` menu item). Search "julia" to find the Julia extension and click Install. A restart of VSCode may then be required. You'll know the Julia extension is working correctly if you can load a `.jl` file and see syntax highlighting. The blue status bar at the bottom of the window will also print some Julia information. The first launch of the Julia extension may take a while to complete (for example, the "Indexing packages..." step might take a couple minutes). Once the extension has fully loaded, a lot of powerful features become available. To see how it should look, see the [Julia for VSCode](https://www.julia-vscode.org/) landing page. Features include "auto-complete" suggestions while typing, pop-up documentation on mouse-hover, an integrated REPL, an integrated debugger, a plot panel, and more.

VSCode does bring a significant learning curve. The Command Palette can help a lot with discoverability. Access it through the `View -> Command Palette...` menu item (Shift-Command-P on Mac, or Shift-Ctrl-P on Windows). Here you can start typing a command to get a list of suggestions. For example, entering "julia repl" will suggest the `Julia: Start REPL` command. Press Enter to launch a Julia REPL, which will appear at the bottom of the VSCode window. This running Julia process is integrated with other Julia-VSCode features in a powerful way.

Again bring up the Command Palette (`View -> Command Palette...`) and type "settings json" to find the command `Preferences: Open Settings (JSON)`. Press enter to and you will an editor for the `Settings.json` file. Add to the bottom of this file the line

```
    "julia.execution.resultType": "inline",
```

Upon saving this settings file, VSCode will immediately adopt the changes.

You can interactively evaluate code in any file `.jl` using the Option-Enter (Mac) or Alt-Enter (Windows/Linux) key command, which maps to `Julia: Execute Code in REPL and Move`. This command will send the Julia expression under the cursor to the running Julia REPL for evaluation. The result of evaluation will be displayed "inline" in the text editor using a distinct color. Effectively, one gets the power and interactivity of a Jupyter notebook, but with the convenience of working with ordinary `.jl` files.

Every window in VSCode represents a standalone process. In most cases, you will probably want to do all work entirely in a single VSCode window. To develop a package, it is useful to load the directory containing all source files into the VSCode window (e.g., using `File -> Open ...`). You can then navigate the files from within VSCode.

It is frequently useful to launch VSCode from the terminal. On Unix systems, run `Shell Command: Install 'code' command in PATH` using the Command Palette. (For Windows systems, no action is necessary.) As expected, this will create a shell command `code` that can be used to quickly launch VSCode. The usage `code <filename>` and `code <directory>` is also supported.

It is convenient to make VSCode the default editor for Julia. On a UNIX system, this is possible by adding the line
```bash
export JULIA_EDITOR=code
```
to the shell startup script (e.g. `.bashrc` or `.bash_profile`). On a Windows system, it [seems the best way](https://discourse.julialang.org/t/windows-command-for-vs-code/29902/9) to configure this environment variable is to add the line `ENV["JULIA_EDITOR"]="code.cmd"` to the `.julia/config/startup.jl` file (note the extra `.cmd` suffix). The file `startup.jl` will not exist in a fresh Julia install; you will need to create it by hand.

Once `JULIA_EDITOR` has been configured, the `@edit` macro can be used to load source code. For example, running
```julia
julia> @edit sort([3, 2, 1])
```
will open the definition of the `sort` function in the VSCode editor. Besides `@edit`, Julia supports several other [interactive utilities](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/).

Another highly recommended VSCode extension is [Git History](https://marketplace.visualstudio.com/items?itemName=donjayamanne.githistory).
