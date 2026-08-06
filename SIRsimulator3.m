% SIRsimulator.m
function [Rnor_all, Rmis_all, Rnor0, Pnor0, Pnor_all, Pmis_all] = SIRsimulator3(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate)
% A function to simulate the spread of misfolded alpha-syn

%% input parameters (inside parenthesis are values used in the paper)
% N_regions: number of regions (42)
% v: speed (1)
% dt: time step (0.01)
% T_total: total time steps (10000)
% GBA: GBA gene expression (zscore, N_regions * 1 vector) (empirical GBA expression)
% SNCA: SNCA gene expression after normalization (zscore, N_regions * 1 vector) (empirical SNCA expression)
% sconnLen: structural connectivity matrix (length) (estimated from HCP data)
% sconnDen: structural connectivity matrix (strength) (estimated from HCP data)
% ROIsize: region sizes (voxel counts)
% seed: seed region of misfolded alpha-syn injection (choose as you like? (^?^)= here substantia nigra)
% syn_control: a parameter to control the number of voxels in which
% alpha-syn may get synthesized (region size, i.e., ROIsize)
% init_number: number of injected misfolded alpha-syn (1)
% prob_stay: the probability of staying in the same region per unit time (0.5)
% trans_rate: a scalar value, controlling the baseline infectivity

%% output parameters
% Rnor_all: A N_regions * T_total matrix, recording the number of normal
% alpha-syn in regions
% Rmis_all: A N_regions * T_total matrix, recording the number of
% misfolded alph-syn in regions
% Pnor_all: a N_regions * N_regions * T_total matrix, recording the number of normal alpha-syn in paths
% could be memory-consuming)
% Pmis_all: a N_regions * N_regions * T_total matrix, recording the number of misfolded alpha-syn in paths
% could be memory-consuming)
% Rnor0: a N_Regions * 1 vector, the population of normal agents in regions before pathogenic spreading
% Pnor0: a N_Regions * 1 vecotr, the population of normal agents in edges before pathogenic spreading

% make sure the diag is zero
sconnDen(1:N_regions+1:end) = 0;
sconnLen(1:N_regions+1:end) = 0;

sconnMov = v ./ sconnLen .* dt;
sconnMov(sconnLen == 0) = 0; % longer path & smaller v = lower probability of moving out of paths

% set the mobility pattern
weights = sconnDen;
weights = (1 - prob_stay) .* weights + prob_stay .* diag(sum(weights, 2)) ;

% multinomial distribution
% element (i,j) is the probability of moving from region i to edge (i,j)
weights = weights ./ sum(weights, 2) .* dt;
weights(1:N_regions+1:end) = 0;

% convert gene expression scores to probabilities
clearance_rate = normcdf(zscore(GBA));
synthesis_rate = normcdf(zscore(SNCA));

clearance = exp(-clearance_rate .* dt);
synthesis = (synthesis_rate .* syn_control).*dt;

gamma0 = trans_rate ./ ROIsize .* dt;

% store the number of normal/misfoled alpha-syn at each time step
[Rnor_all, Rmis_all] = deal( zeros([N_regions, T_total]) );
% [Pnor_all, Pmis_all] = deal( zeros([N_regions, N_regions, T_total]) );

% Rnor, Rmis, Pnor, Pmis store results of single simulation at each time
[Rnor, Rmis] = deal(zeros(N_regions, 1)); % number of normal/misfolded alpha-syn in regions
[Pnor, Pmis] = deal(zeros(N_regions)); % number of normal/misfolded alpha-syn in paths


%% normal alpha-syn growth
% fill the network with normal proteins
N = N_regions;

%%% 1. Build the Markov Transition Components
% To make this efficient and prevent memory explosion, we use sparse matrices.

% [R -> R]: Proteins that stay in regions (1 minus sum of outward weights)
M_RR = spdiags(1 - sum(weights, 2), 0, N, N);

% [P -> P]: Proteins that stay in paths (1 minus outward path connections)
M_PP = spdiags(logical(sconnMov(:)) - sconnMov(:), 0, N^2, N^2);

% For moving between R and P, we map the 2D matrix indices to the 1D vector
[rows_i, cols_j] = ndgrid(1:N, 1:N); 

% [R -> P]: R(i) sends proteins to P(i,j) based on weights(i,j)
% rows_i(:) gives the source region 'i' for every linear index in P
M_RP = sparse(1:N^2, rows_i(:), weights(:), N^2, N);

% [P -> R]: P(i,j) sends proteins to R(j) based on sconnMov(i,j)
% cols_j(:) gives the destination region 'j' for every linear index in P
M_PR = sparse(cols_j(:), 1:N^2, sconnMov(:), N, N^2);

% Combine into a single Movement matrix with clearance in the regions afterwards:
C_diag = spdiags(clearance(:), 0, N, N);
M = [C_diag * M_RR, C_diag * M_PR; 
          M_RP, M_PP];

%%% 2. Synthesis
B = [synthesis(:); zeros(N^2, 1)];

%%% 3. Run the Markov Loop
% Initialize unified state vector X = [Rnor; Pnor(:)]
X = zeros(N + N^2, 1);

iter_max = 1000000000;
display('normal alpha synuclein growth (Markov System)')

for t = 1:iter_max
    X_tmp = X(1:N); % Save old Rnor to check for convergence

    % The entire moving and growth process reduces to one line:
    X = M * X + B;

    if abs(X(1:N) - X_tmp) < (1e-7 * X_tmp)
        break;
    end
end

%%% 4. Unpack results back to original shapes
Rnor0 = X(1:N);
Pnor0 = reshape(X(N+1:end), N, N);

Rnor = Rnor0;
Pnor = Pnor0;


%% misfolded protein spreading process

% inject misfolded alpha-syn
Rmis(seed) = init_number;
display('misfolded alpha synuclein spreading')
for t = 1:T_total
    %%% moving process
    % normal proteins: region -->> paths
    movDrt_nor = Rnor .* weights;

    % normal proteins: paths -->> regions
    movOut_nor = Pnor .* sconnMov;

    % misfolded proteins: region -->> paths
    movDrt_mis = Rmis .* weights;

    % misfolded proteins: paths -->> regions
    movOut_mis = Pmis .* sconnMov;

    % update regions and paths
    Pnor = Pnor - movOut_nor + movDrt_nor;
    Rnor = Rnor + sum(movOut_nor, 1)' - sum(movDrt_nor, 2);

    Pmis = Pmis - movOut_mis + movDrt_mis;
    Rmis = Rmis + sum(movOut_mis, 1)' - sum(movDrt_mis, 2);

    % the probability of getting misfolded
    misProb = 1 - exp( -Rmis .* gamma0) ; % trans_rate: default
    % number of newly infected
    N_misfolded = Rnor .* exp(-clearance_rate) .* misProb ; % is this supposed to be -clearance_rate.*dt?

    % update
    Rnor = Rnor.*clearance - N_misfolded + synthesis;
    Rmis = Rmis.*clearance + N_misfolded;
    Rnor_all(:, t) = Rnor ;
    Rmis_all(:, t) = Rmis ;

    % uncomment the following lines if you want outputs of alpha-syn in
    % paths
    %Pnor_ave(:, :, t) = Pnor;
    %Pmis_ave(:, :, t) = Pmis;

end
end


