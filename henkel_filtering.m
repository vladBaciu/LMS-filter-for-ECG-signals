%% Filtering with reduced Henkel Matrix
% The following method is proposed in the following publication
% Filtering via Rank-Reduced Hankel Matrix, CEE 690 , ME 555 — System Identification — Fall, 2013


N = 2048;
y = ECG_validation_waveform(1:N);
dt = 0.05;
t = [1: N ]* dt ;
y2= y;
m = ceil(0.6* N + 1 );
n = length(y) +1 - m;
Y = zeros (m , n );
sv_ratio = 0.04; 

%% Create a Henkel matrix with the input signal ECG_validation_waveform
for k =1: m
    Y ( k , :) = y ( k :k +n -1 );
end

[U ,S ,V] = svd (Y , 0);

%% Get a maximum index from the matrix of singluar values
f = find(diag(S)/ S(1 ,1) > sv_ratio );
K = max(f);
d1= diag(S)';

%% Create a rank-reduced Hankel Matrix according to K index
d = diag(d1(1:K));
v = V(:,1: K)';
Y = U(:,1: K )*d* v;


y = zeros (1 , N );
y (1) = Y (1 ,1);

%% Filter the signal
for k =2: m
    min_kn = min(k , n );
    first_col = flip(Y(1:1:k ,1: min_kn));
    y(k) = sum(diag(first_col)) / min_kn ;
end


for k =2: n
    last_col = flip(Y(m-n+k:1:m,k:n));
    y(m+k -1) = sum(diag(last_col)) / (n -k +1);
end
figure
plot (t , y2,t,y)