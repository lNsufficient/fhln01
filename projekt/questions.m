% 2017-03-01:
%Hur kan man använda symmetri när den inte är symmetrisk kring
% inspänningar. Svar: Egentligen har vi frihet i x-led vid båda
% inspänningar, men för att förhindra stelkroppsrörelse fixerar vi ett av
% stöden. Dock hade vi i implementeringen inte behövt ha det såhär nu när vi
%har symmetri i mitten. 


% Vad ska vi ha för elasticitetsmodul? %Välj själva, t.ex. för stål

% MMA eller CONLIN? %Vi får välja här - fär gärna även testa MMA, många har
% gjort deta istället för CONLIN. 

%Niklas tycker att vår lösningen ser bra ut och att det funkar att färga in
%elementen som slår i lägre gräns så som vi gjort. Man kan försöka variera
%tolerans och A_min så att den övre högra tirangeln ser bra ut. 
%-------------------------------------------------------------
%Att fundera på: 
%----------------------------------
%Kan vi på något sätt se om våra spänningar ser vettiga ut? 
%Ja - de element som inte slagit i någon gräns ska ha samma spänning, 
%sqrt(lambdastar/E).
    %Varför stämmer då inte detta?

%Vad händer om vi gör en urspriunglig gissning på något som är fel?

%Hur ska man tänka kring tolerans?

%Måste ae uppdateras som följd av deformering?

%(Hur väljs alpha?)

% Vad förväntar man sig för typ av konvergens? På coarse mesh blir det räta
% linjer då semilogy(res) plottas, men på fine mesh så blir de inte räta
% (med positiv andraderivata), så konvergenshastigheten minskas hela tiden.
% Ska det vara så?

%Är det okej att anta att phi_i går att räkna bara på mittpunkten av ett element?