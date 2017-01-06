% VSDP gives the possibility to choose via global parameters
% alternative numerical ways:

% Choice of the SDP-solver in MYSDPS
% VSDP_CHOICE_SDP = 1: SDPT3 solver
% VSDP_CHOICE_SDP = 2: SDPA solver
% DEFAULT = 1;
global VSDP_CHOICE_SDP; VSDP_CHOICE_SDP = 1;

% Choice, whether in VULS the interval system is solved as a full
% system, or as a sparse system by computing the sdp matrix B*B';
% VSDP_CHOICE_FULL = 0 in the latter case, recommended only for very.
% large system. DEFAULT = 1;
global VSDP_CHOICE_FULL; VSDP_CHOICE_FULL = 1;

% Choice, whether during the iteration VSDPLOW and VSDPUP
% the last approximate solution is used as a starting point
% for the next step (VSDP_USE_STARTING_POINT = 1), or not
% (VSDP_USE_STARTING_POINT = 0 is DEFAULT);
global VSDP_USE_STARTING_POINT; VSDP_USE_STARTING_POINT = 0;

% Choice of the maximal number of iterations in VSDP
% DEFAULT = 5;
global VSDP_ITER_MAX; VSDP_ITER_MAX = 9;

% Choice of the perturbation parameter; DEFAULT = 1.5;
global VSDP_ALPHA; VSDP_ALPHA = 1.5;




% Using higher precision for computing the defects in VULS:
% VSDP_HIGHER_PREC = 0 implies usual interval arithmetic, otherwise
% we use the INTLAB function DOT_ with quadruple precision
% DEFAULT = 0;
global VSDP_HIGHER_PREC; VSDP_HIGHER_PREC = 0;
