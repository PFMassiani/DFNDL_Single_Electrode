%% Purpose
% Test the functions defined in this file

%% matrices_linearized_t_indep.m

params = params();

N_elyte = params.dscrtzn.N_e_n;
N_s = params.dscrtzn.N_s_n;

ELYTE = ones(N_elyte+1,1);
ELYTE2N = ELYTE(2:end-1);

US = ones(N_s+1,1);
US2N = US(2:end-1);

USSURF = ones(N_elyte+1,1);
USSURF2N = USSURF(2:end-1);

DLK = 1;


TOLERANCE = 10^(-10);
matrices = matrices_linearized_t_indep(params);

disp 'Testing matrices_linearized_t_indep.m...'
disp '---------------------------------------------------------------'
disp 'Scaled differentiation operators (no BC)'
disp ' '

testvar = matrices.adim.ephase.D1_noBC * ELYTE;
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = 'FAILED';end
msg = ['First order differentiation of constant (electrolyte phase) ....... ', result];
disp(msg)

testvar =  matrices.adim.sphase.D1_noBC * US;
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = 'FAILED';end
msg = ['First order differentiation of constant (solid phase) ....... ' ,result];
disp(msg)

testvar = matrices.adim.ephase.D2_noBC * ELYTE;
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = 'FAILED';end
msg = ['Second order differentiation of constant (electrolyte phase) ....... ', result];
disp(msg)

% !!! A linear function is not linspace distributed, but "chebspace" distributed !
US = chebspace(0,N_s,N_s);
US2N = US(2:N_s);

testvar = matrices.adim.sphase.D2_noBC * US;
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(testvar))];end
msg = ['Second order differentiation of linear function (electrolyte phase) ....... ', result];
disp(msg)

disp ' '
disp 'Boundary conditions-satisfying differentiation operators'
disp ' '

testvar = matrices.adim.ephase.elyte.D2_BC * ELYTE2N + matrices.adim.ephase.elyte.exogeneous_bdry;
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(testvar))];end
msg = ['Second order differentiation of constant (electrolyte phase) ....... ', result];
disp(msg)

USSURF = 3;
US = USSURF * chebspace(0,1,N_s);
US2N = US(2:end-1);
DLK = USSURF;

testvar = matrices.adim.sphase.us.D2_BC * US2N + matrices.adim.sphase.us.exogeneous_linear_dl * DLK - matrices.adim.sphase.us.exogeneous_nonlinear(DLK,US2N)%;
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(abs(testvar)))];end
msg = ['Second order differentiation of linear function (solid phase) ....... ', result];
disp(msg)
disp 'TODO : resolve this issue. Ce qui foire = matrices.adim.sphase.us.exogeneous_linear_dl. Refaire la BC à r = 1'