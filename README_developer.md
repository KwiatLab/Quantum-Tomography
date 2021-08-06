# Developing on Quantum-Tomography
This package is under an MIT License, meaning it is free to be expanded upon as long as it is kept open source. 
For more details see our LICENSE file.

### Contents
- [Setting up your Environment](#setting-up-your-environment)
- [Branching Strategy](#branching-strategy)
- [Publishing a new Version](#publishing-a-new-version)
- [Updating Documentation](#setting-up-your-environment)
- [Github Actions](#github-actions)



## Setting up your environment
This python application is tested on many python version. It is recommended to develop on 
the [most recent version of python](https://www.python.org/downloads/). Tests are setup to run on previous versions.

1. **Clone** the repo to your local computer


2. **Create a virtual environment** by typing the following in the command prompt. Make sure you are at the top most level of the repo.

         py -m venv src\venv
3. **Activate virtual environment**

        src\venv\Scripts\activate.bat

4. Install the **requirements**.

         pip install -r requirements.txt
         
5. **Switch to the development branch** before making any changes
         
        git checkout development

If you are using pycharm(a great development environment for python), you can add the virtual environment by going to File->Settings->Project->python interpreter. 
   Click the drop down for interpreter and click show all. Click the + button, and add a new existing interpreter using the path:
   _C:\path\to\local\repo\src\venv\Scripts\python.exe_. Any IDE can be used to develop this code.


### Troubleshooting

- *'py' is not recognized*. 
  - Depending on how you downloaded python your PATH variable for python (the _py_ in this walkthrough) may be diffrent.
Other common names are _python_ or _python3_. Its whatever you type into the command prompt in order to start python. 
  Replace the _py_ with whatever your PATH variable for python is. Follow [this tutorial](https://www.educative.io/edpresso/how-to-add-python-to-path-variable-in-window) 
  if you can't figure out your python's PATH variable and you know for certain python is installed.

- *'pip' is not recognized*. 
  - pip is not defined as a PATH variable. Instead of just using *pip*
use *py -m pip*. It is rare but if you downloaded python in a weird way you may not have pip installed. You'll have to search for how to install pip on your computer. 
The process differs depending on if you're using PC/mac/ubuntu.
    

## Branching Strategy
We have adopted **Git-Flow** as our branching strategy. This is a set of rules
to develop on the code. Git-Flow is a very popular strategy used in the field. Here
is a [video going over Git-Flow](https://youtu.be/y4yg7aT4NgM?t=80), however heres a quick summary:

#### master
This code is always depolyable. In theory we could setup our github actions script to publish our
code everytime we push to master. However instead of using a **Release** branch, we use master branch in case there
are any last minute issues. When it is ready to go, we push the current state of master to pypi.

#### development
This is the main branch you will be making changes to. When ready you can create a pull request from development into master.

#### feature
For some changes you may want to branch off of development, make edits, then merge back into development. 
You want to avoid your feature branch from diverging from development. Here are some tips:
 - Keep your feature branches short and merge fast.
 - When possible and sensible, make changes to the development then merge those changes into your feature branch.
    - ex: Alice and Bob both branch off of development and are creating new functions f_A and f_B. Bob wants to rename a variable in 
      the filter_data function. Bob should make this change in the development branch then both Bob and Alice should pull that change into
      their respected featuer branches.
  
#### hotfix
If you need to make a small change to the production code, but there are features in the development branch that are not 
ready to go, use a hotfix. Branch off of master, make the change, and merge this into BOTH master and development (and any feature branches)

## Publishing a new Version
Please see our Group's Wiki if you are looking to publish a new version.

## Updating Documentation
Documentation is automatically created by our flask site. Comment blocks above functions are parsed into html pages.
The format of the comment blocks is provided in the __init__.py script in the src folder. TIP: If you want to exclude a function's comment block
from being displayed on the flask site, use single quotes ''' instead of double quotes """.


# Github Actions
We have two github actions setup. One is used to publish the code to Pypi and the other is used to run
our test scripts across multiple versions. These actions are defined in the .github/workflow folder. The associated tests
that are run with each actions are in the Tests folder titled test_github_*.py. SEE ALSO : README.md in Tests.

