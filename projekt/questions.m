% 2017-03-01:
%Hur kan man anv�nda symmetri n�r den inte �r symmetrisk kring
% insp�nningar. Svar: Egentligen har vi frihet i x-led vid b�da
% insp�nningar, men f�r att f�rhindra stelkroppsr�relse fixerar vi ett av
% st�den. Dock hade vi i implementeringen inte beh�vt ha det s�h�r nu n�r vi
%har symmetri i mitten. 


% Vad ska vi ha f�r elasticitetsmodul? %V�lj sj�lva, t.ex. f�r st�l

% MMA eller CONLIN? %Vi f�r v�lja h�r - f�r g�rna �ven testa MMA, m�nga har
% gjort deta ist�llet f�r CONLIN. 

%Niklas tycker att v�r l�sningen ser bra ut och att det funkar att f�rga in
%elementen som sl�r i l�gre gr�ns s� som vi gjort. Man kan f�rs�ka variera
%tolerans och A_min s� att den �vre h�gra tirangeln ser bra ut. 
%-------------------------------------------------------------
%Att fundera p�: 
%----------------------------------
%Kan vi p� n�got s�tt se om v�ra sp�nningar ser vettiga ut? 
%Ja - de element som inte slagit i n�gon gr�ns ska ha samma sp�nning, 
%sqrt(lambdastar/E).
    %Varf�r st�mmer d� inte detta?

%Vad h�nder om vi g�r en urspriunglig gissning p� n�got som �r fel?

%Hur ska man t�nka kring tolerans?

%M�ste ae uppdateras som f�ljd av deformering?

%(Hur v�ljs alpha?)

% Vad f�rv�ntar man sig f�r typ av konvergens? P� coarse mesh blir det r�ta
% linjer d� semilogy(res) plottas, men p� fine mesh s� blir de inte r�ta
% (med positiv andraderivata), s� konvergenshastigheten minskas hela tiden.
% Ska det vara s�?

%�r det okej att anta att phi_i g�r att r�kna bara p� mittpunkten av ett element?