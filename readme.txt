This project is a satellite of <supy> for SusyCAF format TTrees.
See <github.com/elaird/supy>.

-----------
| License |
-----------
GPLv3 (http://www.gnu.org/licenses/gpl.html)

---------------
| Quick Start |
---------------
If you have access to CERN AFS, you can try the example.
<set up pyROOT>                               #1) eg. with a CMSSW area: cd /<somepath>/CMSSW_4_2_8/src && cmsenv
git clone git://github.com/<user>/susycaf.git #2) clone repo: <user> can be betchart or a forking user
cd susycaf
git submodule update --init                   #3) checkout supy dependence
export PYTHONPATH=$PYTHONPATH:`pwd`           #4a) add directory containing supy to your python path
export PATH=$PATH:`pwd`/supy/bin	      #4b) optionally add to your path
supy analyses/example.py --loop 1             #5) Run the example (the example input files are located on AFS):

----------------
| Dependencies |
----------------
ROOT (>=5.27.06) and python (2.x, x>=6) are required; CMSSW is not.
See <github.com/elaird/supy>

--------
| Bugs |
--------
Please report problems at https://github.com/betchart/susycaf/issues
