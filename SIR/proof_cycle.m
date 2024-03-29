% This script proves Theorem 12. It will request some inputs
% related to parallel processing. You must have INTLAB installed or the
% script will fail. Note that this does NOT require MATLAB's parallel
% processing toolbox, and in fact, INTLAB is not compatible with
% it anway. You will need to manually run multiple MATLAB instances. 
% Doing it this way is thread safe.
clear all
load('candidate_cheb_N16.mat')
parallel(2) = input('How many threads (n) do you want to run? For a balanced session, choose n that divides 10000. Enter a positive integer. n = ');
parallel(1) = input('Which thread (m) is this?. Enter an integer 1<=m<=n. m = ');
[p,radii] = verify_branch(SI,times,alpha,para,N,M,nu,parallel,W,rstar);
