# Getting started with Julia for users of Sunny

## Julia basics

The latest stable release of [Julia](https://julialang.org/) is available from the [downloads page](https://julialang.org/downloads/). Alternatively, the tool [`juliaup`](https://github.com/JuliaLang/juliaup) offers a convenient way to update Julia when new versions are released.

Once installed, execute the `julia` command, bringing a prompt like this:

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.8.0 (2022-08-17)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
```

This is the Julia REPL, short for read-eval-print-loop. You can run code here, for example
```julia
julia> exp(π * im)
-1.0 - 1.2246467991473532e-16im
```
which verifies Euler's identity `e^{π i} = -1` to floating point accuracy, `≈ 10^-16`.

There are many resources to learn Julia. Quick tutorials are given in [_Learn X in Y minutes, where X=Julia_](https://learnxinyminutes.com/docs/julia/) and the [Julia Cheat-sheet](https://juliadocs.github.io/Julia-Cheat-Sheet/). There is also the [official Julia documentation](https://docs.julialang.org/), but you don't need to be a Julia expert get started with Sunny.


## The built-in Julia package manager

You may be familiar with package managers like `pip` and `anaconda` for Python. In Julia this functionality is built-in.

At the Julia REPL, press the key `]` to enter the "package" mode:

```
(@v1.8) pkg> 
```

New packages can be added with the `add` command. To install Sunny, use:

```
pkg> add Sunny
```

It will be useful to install two plotting packages:
```
pkg> add Plots, GLMakie
```

Many other packages are available, and can be browsed at [JuliaHub](https://juliahub.com/ui/Packages).

Type `help` to see a list of available commands in package mode. For example, to see the list of currently installed packages, enter `status`. Packages can be removed with the `remove` command, or updated with the `update` command.

The [Backspace] key exits package mode.

## Getting help

Julia has built-in help. From the REPL, enter `?` to enter the help mode, and then enter a function or type name to get documentation. For example:
```
help?> exp
search: exp exp2 Expr expm1 exp10 export exponent expanduser ExponentialBackOff ldexp frexp nextpow

  exp(x)

  Compute the natural base exponential of x, in other words e^x.
  ...
```

Public functions of Sunny should be similarly documented.

Math symbols can be accessed with Latex-like notation. For example, typing `\pi` and then pressing Tab will produce `π`. If you see a unicode symbol that you don't know how to type, you can copy-paste it into the help prompt:

```
help?> ≈
"≈" can be typed by \approx<tab>
...
```

Although Julia and Matlab have similar surface syntax, there are important differences between the two languages. One crucial thing to know is that simple "for loops" over scalars are _encouraged_ in Julia, and can be used to generate code with C++ or Fortran-like speed. Vectorized style is supported in Julia, but not necessary. Like Matlab and Fortran (but unlike Python or C++), Julia uses 1-based indexing by default, so `my_array[1]` is the first element of `my_array`.


## Interactive "notebook" environment

We now recommend the [Julia VSCode extension](https://www.julia-vscode.org/) for interacting with Sunny. Julia VSCode enables a notebook-like experience overlayed on top of a full-featured code-editor. Tips for working with Julia VSCode are provided in a later section.

Jupyter notebooks are another option. To enable this integration, one can install the IJulia package,

```
pkg> add IJulia
```

Launch a Julia/Jupyter notebook using:
```julia
using IJulia
notebook()
```

The Julia/Jupyter integration seems to break whenever Julia is updated. If you see an error about a 'missing kernel', reinstall the Julia kernel files by executing `pkg> build IJulia`.


# Advanced features for code development

For beginning users, the previous sections are fully sufficient for running and modifying the Sunny examples. Some users may want to experiment with changes to Sunny itself. This section describes a more advanced setup for Julia code development.

## Developing the Sunny source code

The [Git version control system](https://git-scm.com/) makes it possible to fearlessly experiment with code modifications. To get the latest development branch of Sunny, use:

```
pkg> rm Sunny
pkg> dev Sunny
```

This will download (more specifically, `git clone`) the source code for Sunny into the local directory `~/.julia/dev/Sunny/`. Modifications to the files here will be picked up by Julia.

A full introduction to Git is beyond the scope of this document, but here are some basics. Open a terminal in the `~/.julia/dev/Sunny/` directory and type `git status`. You should see that the directory is free of changes, i.e., "clean". Try modifying some source file. Now `git status` will show name of the file that was changed. Type `git diff` to see the specific changes made. You can revert these changes with the command `git checkout <filename>`. Other useful commands include `git add <filenames>` and `git commit`, which will enter changes into the database (repository) of tracked changes (commits). `git log` will show a history of commits. The commands `git pull` and `git push` will download and upload, respectively, from the "origin" repository (in this case, the one hosted on Github). If you actually try this, you will likely find that `git push` reports an error stating that you don't have write access to the main Sunny repository. The recommend workflow is to [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) Sunny on Github. The forked version can be checked out using `pkg> dev <GithubURL>`. Changes in a fork can be considered for incorporation into Sunny using a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/).

Instead of entering Git commands in the terminal, it's usually more convenient to use a graphical user interface. For Julia development, we recommend the VSCode editor with the Git Lens extension (see below).

## Julia development with Revise.jl

An extremely important package in the Julia ecosystem is [Revise](https://timholy.github.io/Revise.jl/stable/):
```
pkg> add Revise
```

Revise allows a running Julia process to dynamically pick up changes to package source code (i.e., update function behaviors) without requiring a restart. This is very convenient for rapid iteration -- one can avoid most of the wait times associated with reloading (and recompiling) packages.

Enable Revise by default in every Julia session by [following the instructions here](https://timholy.github.io/Revise.jl/stable/config/).

Revise works great for redefining function definitions, but has [some limitations](https://timholy.github.io/Revise.jl/stable/limitations/). In particular, Revise does _not_ support redefining `struct` datatypes. If a `struct` is modified, then the REPL must be restarted.

## The Julia extension in VSCode

[VSCode](https://code.visualstudio.com/) is a powerful text editor, and hosts the official [Julia development environment](https://www.julia-vscode.org/).

Much of the power of VSCode comes from extensions. Download these through the Extensions panel, accessible by clicking the appropriate icon on the left of the VSCode window (alternatively, by selecting the `View -> Extensions` menu item). Search "julia" to find the Julia extension and click Install. A restart of VSCode may then be required. You'll know the Julia extension is working correctly if you can load a `.jl` file and see syntax highlighting. The blue status bar at the bottom of the window will also print some Julia information. The first launch of the Julia extension may take a while to complete (for example, the "Indexing packages..." step might take a couple minutes). Once the extension has fully loaded, a lot of powerful features become available. To see how it should look, see the [Julia for VSCode](https://www.julia-vscode.org/) landing page. Features include "auto-complete" suggestions while typing, pop-up documentation on mouse-hover, an integrated REPL, an integrated debugger, a plot panel, and more.

The Command Palette can help to discover VSCode functionality. Access it through the `View -> Command Palette...` menu item (Shift-Command-P on Mac, or Shift-Ctrl-P on Windows). Here you can type keywords to get a list of command suggestions. For example, entering "julia repl" will suggest the `Julia: Start REPL` command. Press Enter to launch a Julia REPL, which will appear at the bottom of the VSCode window. This running Julia process integrates with other Julia-VSCode features in a powerful way.

Again bring up the Command Palette and type "settings json" to find the command `Preferences: Open Settings (JSON)`. Press enter to and you will an editor for the `Settings.json` file. Add to the bottom of this file the line

```
    "julia.execution.resultType": "inline",
```

Upon saving this settings file, VSCode will immediately adopt the changes.

You can interactively evaluate code in any file `.jl` using the Shift-Enter (Mac), which maps to `Julia: Execute Code in REPL and Move`. This command will send the Julia expression under the cursor to the running Julia REPL for evaluation. The result of evaluation will be displayed "inline" in the text editor using a distinct color. Effectively, one gets the power and interactivity of a Jupyter notebook, but with the convenience of working with ordinary `.jl` files.

Every window in VSCode represents a standalone process. In most cases, you will probably want to do all work entirely in a single VSCode window. To develop a package, it is useful to load the directory containing all source files into the VSCode window (e.g., using `File -> Open ...`). You can then navigate the files from within VSCode.

It is frequently useful to launch VSCode from the terminal. On Unix systems, run `Shell Command: Install 'code' command in PATH` using the Command Palette. On Windows systems, `code` will be available by default. As expected, this will create a shell command `code` that can be used to quickly launch VSCode. The usage `code <filename>` and `code <directory>` is also supported.

It is convenient to make VSCode the default editor for Julia. On a UNIX system, this is possible by adding the line
```bash
export JULIA_EDITOR=code
```
to the shell startup script (e.g. `.bashrc` or similar). On a Windows system, it [seems the best way](https://discourse.julialang.org/t/windows-command-for-vs-code/29902/9) to configure this environment variable is to add the line `ENV["JULIA_EDITOR"]="code.cmd"` to the `.julia/config/startup.jl` file (note the extra `.cmd` suffix). The file `startup.jl` will not exist in a fresh Julia install; you can create it by hand.

Once `JULIA_EDITOR` has been configured, the `@edit` macro can be used to load source code. For example, running
```julia
julia> @edit sort([3, 2, 1])
```
will open the definition of the `sort` function in the VSCode editor. The `@edit` macro is defined in the [`InteractiveUtils` module](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/). Browse around these docs to see what else is availabe.
