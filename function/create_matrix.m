%Corrected for lower weights between EI    
function mat = create_matrix(E,I,EE,EI,IE,II,R,gamma)
    tot = E+I;
    EEmat = create_submatrix(E, E, tot, EE, R, gamma, 0);
    EImat = create_submatrixEI(E, I, tot, EI, R, gamma, 1);
    IEmat = create_submatrix(I, E, tot, IE, R, gamma, 1);
    IImat = create_submatrix(I, I, tot, II, R, gamma, 1);

    Emat = [EEmat;EImat];
    Imat = [IEmat;IImat];
    mat = [Emat Imat];

end

function submat = create_submatrix(N, M, tot, p, R, gamma , inh)

NM = round(p*N*M);
submat = [ones(1,NM), zeros(1,N*M - NM)];
submat = reshape(submat(randperm(N*M)),M,N); 

w0 = sqrt(2)*R/(sqrt(p*(1-p)*(1 + gamma^2)));
submat = submat*(w0/sqrt(tot));

if inh == 1
    submat = -gamma*submat;
end

submat = triu(submat,1) + tril(submat,-1);

end

function submat = create_submatrixEI(N, M, tot, p, R, gamma , inh)

NM = round(p*N*M);
submat = [ones(1,NM), zeros(1,N*M - NM)];
submat = reshape(submat(randperm(N*M)),M,N); 

w0 = sqrt(2)*R/(sqrt(p*(1-p)*(1 + gamma^2)));
submat = submat*(w0/sqrt(tot));

if inh == 1
    submat = gamma*submat*0.5;
end

submat = triu(submat,1) + tril(submat,-1);

end