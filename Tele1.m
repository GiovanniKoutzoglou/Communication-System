%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Εργασία στα Τηλεπικοινωνιακά Συστήματα 3                                      %
% Περιγραφή :                                                                   %    
%   Έχουμε ένα τηλεπικοινωνιακό σύστημα, το οποίο απαρτίζεται από ένα πομπό,    %
%   ένα κανάλι με απόκριση συχνότητας :                                         %
%               { 1 , if |f +/- fo| <= W                                        %
%      C(f) =   {                               ,                               %      
%               { 0 , else                                                      %       
%   προσθήκη AWGN (n(t), ΦΠΙ: No/2.), ένα δέκτη και κατόπιν ένα δειγματολήπτη.  %      
%   Για τη µετάδοση χρησιµοποιούµε διαµόρφωση M−PAM, δηλαδή µεταδίδουµε         %      
%   σύµβολα της µορφής am = ±1.                                                 %        
%   Ο σκοπός της εργασίας είναι να µελετηθεί η επίδοση του συστήµατος, όταν     %  
%   αυτό σχεδιάζεται µε τις παρακάτω επιλογές :                                 %
%   1. Επιλογή των ϕίλτρων ποµπού και δέκτη έτσι ώστε στον δέκτη να έχουµε      %
%       µηδενική διασυµβολική παρεµβολή.                                        %
%   2. Επιλογή των ϕίλτρων ποµπού και δέκτη έτσι ώστε στον δέκτη να έχουµε      %
%       ελεγχόµενη διασυµβολική παρεµβολή µε διπλοδυαδικό παλµό, και σύµβολο    %
%       προς σύµβολο ϕώραση µε προκωδικοποίηση στον ποµπό.                      %
%   3. Επιλογή των ϕίλτρων ποµπού και δέκτη έτσι ώστε στον δέκτη να έχουµε      % 
%      ελεγχόµενη διασυµβολική παρεµβολη µε διπλοδυαδικό παλµό, και εφαρµογή    % 
%      του αλγορίθµου Viterbi στη ϕώραση.                                       %                                        
%                                                                               %
%    Ioannis. M. Koutzoglou, Dec. 2021                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;


%% -------------------- First part of the Project ! -------------------- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         "System Properties"                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 10000;         %% Number of bits
M = 2;             %% Modulation order (2 for 2-PAM)
k = log2(M);       %% Number of bits per symbol
sps = 10;          %% Samples per symbol
EbNo = 25;
alpha = 1;
snr = EbNo + 10*log10(k) - 10*log10(sps);
shape = 'sqrt';
fs = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         "Bits Generation"                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bits = randi([0 1],n,1);
fprintf('bits: %d \n', bits);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         "M-PAM Modulation"                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symbols = pammod(bits,M);
fprintf('\n');
fprintf('symbols: %d \n', symbols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         "AWGN Channel"                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hn = awgn(symbols,snr);
fprintf('\n');
fprintf('hn = %f \n', hn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         "Noise n(t)"                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt = hn - symbols;
fprintf('\n');
fprintf('nt = %f \n', nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          "Gr Filter"                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defining the sinc filter
sincNum = sin(pi*(-fs:1/fs:fs)); % numerator of the sinc function
sincDen = (pi*(-fs:1/fs:fs)); % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

cosNum = cos(alpha*pi*(-fs:1/fs:fs));
cosDen = (1-(2*alpha*(-fs:1/fs:fs)).^2);
cosDenZero = abs(cosDen)<10^-10;
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;
Gr = sincOp.*cosOp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   "Conv of noise and filter"                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr = conv(nt,Gr,'same');
fprintf('\n');
fprintf('nr = %f \n', nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      "Received symbols"                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yn = symbols + nr;
receiveSym = zeros(size(yn));
for i = 1:length(yn)
      if yn(i) > 0
          receiveSym(i) = 1;
      else
          receiveSym(i) = -1;
      end
end
fprintf('\n');
fprintf('Receive symbols: %d \n', receiveSym);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          "Error Rate"                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error = symbols ~= receiveSym;
sum_error = sum(error);
error_prob = (sum_error / n)*100;
fprintf('\n');
fprintf('error prob : %f \n',error_prob);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          "Plots"                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Symbols plot
figure;
stem(real(symbols),':','markerfacecolor', 'black');
title('Symbols');
grid on;
set(gca,'ylim',[-1.5 1.5])

% AWGN Channel plot
figure;
stem(real(hn),':','markerfacecolor','green')
title('hn');
grid on;
set(gca,'ylim',[-1.5 1.5])

% Noise nt plot
figure;
stem(real(nt),':','markerfacecolor','black')
title('nt');
grid on;
set(gca,'ylim',[-1.5 1.5])

% Gr Filter plot
figure;
stem((-fs:1/fs:fs),(Gr),'b','LineWidth',2)
grid on
xlabel('time, t')
ylabel('amplitude, g(t)')
title('Time domain waveform of raised cosine pulse shaping filters')

% Conv of noise and filter plot
figure;
stem(real(nr),':','markerfacecolor', 'magenta');
title('nr(t)');
grid on;
set(gca,'ylim',[-1.5 1.5])

% Receiver plot
figure;
stem(receiveSym,':','markerfacecolor','black')
title('Received Symbols');
grid on;
set(gca,'ylim',[-1.5 1.5])

