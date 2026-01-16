# hello-world

**If you've setup your development environment go to [Assignment](#assignment)**

## Requirements

- Git
- C++ Compiler (gcc, clang, msvc)
- CMake

## Installation

Setting up a C++ development environment can be non-trivial. The course is set up so that most relatively standard development environments will work fine, whether you are using Windows or Linux. Nor do we prescribe that you use a specific editor or that you use an editor at all! However for simplicity here we give instructions on how to use VS Code as your editor and setup a C++ development in Linux or Windows. Microsoft provides a thorough [set of instructions](https://code.visualstudio.com/docs/languages/cpp), but here we will go through some basics and the most pertinent set of instructions.

### Terminal Basics

Many of the following instructions will require use of a terminal. While it is possible to perform these operations purely within a Graphical User Interface (GUI), you will find that you will often run into programs and tools that have no UI, and must be run via a command line interface (CLI), or even if there is a UI, it would be faster to type the necessary commands. The following instructions will be for Windows, as the specific instructions can vary greatly depending on your Linux Distro.

If you are unused to using the terminal here is a brief [introduction](https://www.freecodecamp.org/news/command-line-for-beginners/).

### Setup

#### Windows

On Windows for this first part you can use either Command Prompt or Powershell as your terminal, which can both be accessed via `Windows Terminal`. Both are acceptable but I would suggest using Powershell.

We will use the utility `winget` to install the required programs.

```shell
winget --version
```

First, we will install [Git](https://git-scm.com/), a version control system.

```shell
winget install Git.Git
```

Next, we will install a C++ development environment. On Windows, developers have the option to either use [MSYS2](https://www.msys2.org/), which creates a Unix-like development environment, [Windows-Subsystem-for-Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install), which creates an actual Linux system that runs on top of Windows, or [Microsoft's C++ Build Tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022). If you wish to develop as if you are on Linux, WSL is an option, but we will continue with the Build Tools. While all three can be used with VS Code, the build tools can also be used with Microsoft's [Visual Studio](https://visualstudio.microsoft.com/downloads), which provides a fully integrated development environment. While we're at it, we will also install [VS Code](https://visualstudio.microsoft.com/downloads)

```shell
winget install Microsoft.VisualStudio.2022.BuildTools
winget install Microsoft.VisualStudioCode
```

With the installation of the build tools, you now have two new terminals installed. These are the `Developer Powershell for VS 2022` and `Developer Command Prompt for VS 2022`. The C++ development tools are only available while using one of these new terminals. Open one of the terminals and check that your tools are properly installed

```shell
cl
```

```shell
code -v
```

With everything setup we will make a folder to place some new projects, you could do this in Windows Explorer, but let's continue using the terminal. We will also get the this repository

```shell
cd ~
mkdir Repos/School/mctr294
cd Repos/School/mctr294
git clone https://github.com/mctr294-2026/lab-1-UserName helloWorld
code helloWorld
```

#### Linux

On Linux, use your distribution's package manager to install the required tools. The examples below use `apt` (for Debian/Ubuntu-based systems), but you can substitute with `dnf` (Fedora), `pacman` (Arch), or your system's equivalent.

First, update your package list:

```shell
sudo apt update
```

Install Git, g++ compiler, CMake, and VS Code:

```shell
sudo apt install git g++ cmake
```

For VS Code, you can either download it from the [official website](https://code.visualstudio.com/) or install it via snap:

```shell
sudo snap install --classic code
```

Verify the installations:

```shell
git --version
g++ --version
cmake --version
code --version
```

With everything setup, create a folder for your projects and clone this repository:

```shell
cd ~
mkdir -p Repos/School/mctr294
cd Repos/School/mctr294
git clone https://github.com/mctr294-2026/lab-1-UserName helloWorld
code helloWorld
```

## Building

```shell
cmake -S . -B build
cmake --build build --config Debug
build\Debug\hello_world.exe
```

## Testing

[Build](#building)

```shell
ctest --test-dir build -C Debug
Add -V for verbose testing
```

## Assignment

After you have finished these steps switch to the `roots` branch for part 2 of this lab

``` shell
git switch roots
```

