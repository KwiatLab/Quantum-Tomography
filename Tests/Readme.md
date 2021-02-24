# Testing Foulder
This folder is used to run tests on the src code to ensure
it is working right. It also contains scripts to run random tomography and save results

#  Test on the src code


## Testing
Test.py is a script that runs random tomographies with specified settings.

#### Github
test_github_publish.py and test_github_update.py are used by github
actions. They are automatically run when a new version of a code is published or on every a new push(update). 
They both use the Test.py function

## Saving
SaveRun.py is a class that runs tomographies and saves the results and saves the results. It can generate html with 
all the results of the tomographies. It can also print out a graph of the fidelities. The save_RandomStates.py uses this class.
#### Results
Results are not tracked by git. However the css file in the results folder is tracked. This is used to make the display
of the tests look better in html. 
#### Test_States
This a foulder of data to run through the tomography. Test cases with predetermined counts. It does
not use the SaveRun class.