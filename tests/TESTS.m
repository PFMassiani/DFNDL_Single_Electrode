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

disp 'Now testing matrices_linearized_t_indep.m...'
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

testvar = matrices.adim.ephase.elyte.D2_BC * ELYTE2N + matrices.adim.ephase.elyte.exogeneous_bdry(0);
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(testvar))];end
msg = ['Second order differentiation of constant (electrolyte phase) ....... ', result];
disp(msg)

US = 3 * chebspace(0,1,N_s);
US2N = US(2:end-1);
USSURF = US(1);
DLK = USSURF;

testvar = matrices.adim.sphase.us.D2_BC * US2N + matrices.adim.sphase.us.exogeneous_nonlinear(USSURF);
if (all(testvar < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(abs(testvar)))];end
msg = ['Second order differentiation of linear function (solid phase) ....... ', result];
disp(msg)

disp '---------------------------------------------------------------'
disp 'Test of matrices_linearized_t_indep.m : FINISHED'
disp 'Now testing matrices_linearized_t_dep.m...'
disp '---------------------------------------------------------------'
disp ' '

DL = matrices.adim.ephase.x .* (matrices.adim.ephase.x - 1);

testvar = matrices.adim.ephase.D2_noBC * DL;
if (all(abs(testvar-2) < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(abs(testvar-1)))];end
msg = ['Second order differentiation of parabola (no BC) (electrolyte phase) ....... ', result];
disp(msg)

T = 1:3;
matrices = matrices_linearized_t_dep(matrices,params);

for t = T
    %Building a function DL that satisfies the BC : a parabola, whose
    %derivative in 0 and 1 is equal to -/+ i.
    i = params.adim.i(t);
    CONST = t; % Arbitrary constant : should not change the results
    DL = matrices.adim.ephase.x .* (matrices.adim.ephase.x - 1) * i + CONST;
    DL2N = DL(2:end - 1);
    testvar = matrices.adim.ephase.dl.D2_BC(t) * DL2N + matrices.adim.ephase.dl.exogeneous_bdry(t);
    if (all(abs(testvar - 2*i) < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(abs(testvar - 2*i)))];end
    msg = ['Second order differentiation of parabola (electrolyte phase), TIMESTEP = ',num2str(t),' ....... ', result];
    disp(msg)
end

x = chebspace(0,1,N_elyte);
x = x(2:end-1);
ELYTE = (x.^3)/3 - (x.^2)/2;

testvar = matrices.adim.ephase.elyte.D2_BC * ELYTE + matrices.adim.ephase.elyte.exogeneous_bdry(0);
if (all(abs(testvar-(2*x-1)) < TOLERANCE));result = 'PASSED' ; else; result = ['FAILED, max value = ',num2str(max(abs(testvar-1)))];end
msg = ['Second order differentiation of third order polynomial (electrolyte phase) ....... ', result];
disp(msg)