# Testing Folder

This folder is used to run tests on the src code to ensure
it is working right. It also contains scripts to run random tomography and save results

#  Test on the src code
You should create a virtual environment. It's different depending on your operating system but the command
should be something similar to "py -m venv env"  Once you've done this you'll want to install the local code
as a python package. To do this is the command "py -m pip install -e .". Now when you run the scripts in Tests they
will be referencing the src code in your local directory. If this is not done it's possible the
tests are running on the published code that is downloaded over pypi.

## Testing
runTests()in the TestRun.py is a function that runs random tomographies with specified settings. This returns
a list of tomography objects as well as some other details about how each tomography went.

#### Github
test_github_publish.py and test_github_update.py are used by github
actions. They are automatically run when a new version of a code is published or on every a new push(update).

#### Results
Results files are not tracked by git. These are csv files containing data on the tomographies. The analyze_*.rmd
are R markdown files used to analyze the results and output a html file to view them.

#### Test_States
This a foulder of data to run through the tomography. Test cases with predetermined counts.