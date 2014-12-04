% function do_PC_level4 extracts Principal Components at various levels of the hierarchy

function [F,G,H,FH,FX,GX,K_F,K_G,K_H]=do_PC_level4(bigX,Bsub,k_max,K_F,K_G,K_H);

B = length(Bsub);
bigXmat = []; bigGmat = []; bigHmat = [];

for b=1:B,
    Xb = bigX{b};
    if iscell(Xb), 
        Xb=cell2mat(Xb);
    end;
    bigXmat=[bigXmat Xb];
    [ic_Gb,chat_Gb,GX_b,eigval_Gb]=nbplog(Xb,k_max,2,2);
    if isempty(K_G), K_G(b) = min(max(ic_Gb,1),3); end;
    [ehat_Gb,GX{b},lambda_Gb,ss_Gb]=pc_T(standard(Xb),K_G(b));
    for s = 1:Bsub(b),
        Zbs = bigX{b}{s};
        [ic_Hbs,chat_Hbs,HX_bs,eigval_Hbs]=nbplog(Zbs,min(size(Zbs,2),k_max),2,2);
        if isempty(K_H), K_H{b}(s) = min(max(ic_Hbs,1),3); end;
        [ehat_Hbs,H{b}{s},lambda_Hbs,ss_Hbs]=pc_T(standard(Zbs),K_H{b}(s));
        bigHmat=[bigHmat H{b}{s}];
    end;
    if Bsub(b) > 0,
        [ic_Gb,chat_Gb,GX_b,eigval_Gb]=nbplog(Xb,k_max,2,2);
        if isempty(K_G), K_G(b) = min(max(ic_Gb,1),3); end;
        [ehat1,G{b},lambda_Hbs,ss2]=pc_T(standard(cell2mat(H{b})),K_G(b));
    else,
        G{b} = GX{b};
        K_H{b} = 0;
    end;
    bigGmat=[bigGmat G{b}];
end; % end B

% Extract F from H's
[ic_FH,chat_FH,fhat_FH,eigval_FH]=nbplog(bigHmat,min(k_max,size(bigHmat,2)),2,2);
[ehat2,FH,lambda_FH,ss3]=pc_T(bigHmat,ic_FH);

% Extract F from X's
[ic_FX,chat_FX,fhat_FX,eigval_FX]=nbplog(bigXmat,k_max,2,2);
if isempty(K_F), K_F = min(max(ic_FX,k_max),3); end;
[ehat3,FX,lambda_FX,ss5]=pc_T(standard(bigXmat),K_F);

% Extract F from G's
[ic_F,chat_F,fhat_F,eigval_F]=nbplog(bigGmat,min(k_max,size(bigGmat,2)),2,2);
[ehat2,F,lambda_F,ss4]=pc_T(bigGmat,min(max(ic_F,1),K_F));



