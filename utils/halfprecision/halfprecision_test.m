% halfprecision_test is a test function for halfprecision function
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    halfprecision_test
%  Filename:    halfprecision_test.m
%  Programmer:  James Tursa
%  Version:     2.0
%  Date:        May 18, 2020
%  Copyright:   (c) 2020 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
%******************************************************************************

function halfprecision_test
disp('-----------------------------------------------------------------------');
disp('-----------------------------------------------------------------------');
disp(' ');
disp('                halfprecision( ) test function');
disp(' ');
disp('This function runs different tests depending on the MATLAB version:');
disp(' ');
disp('R2019a and earlier: Runs limited spot checks for expected results');
disp('R2019b and later:   Runs extensive tests compared to MATLAB half( )');
disp(' ');
fprintf('Version:  %s\n',version);
disp(' ');
try
    v = hex2dec(version('-release'));
catch
    v = 0;
end
halfprecision(0); % Autocompile it necessary
halfprecision(0); % Load into memory
R2019b = hex2dec('2019b');
if( v < R2019b )
    halfprecision_test_2019a;
else
    halfprecision_test_2019b;
end
R2017b = hex2dec('2017b');
if( v > R2017b && hex2dec(halfprecision('version')) == R2017b )
   disp(' ');
   w = {'halfprecision( ) has been compiled with R2017b memory model\n',
        'but it is being used in R2018a or later. This will cause\n',
        'deep data copies of complex variables to occur. While this\n',
        'will work, it would be best to recompile halfprecision.c using\n',
        'the -R2018a mex option to avoid these unnecessary deep copies.%s\n'};
   warning([w{:}],' ');
   disp(' ');
end
if( halfprecision('openmp') == 0 )
   disp(' ');
   disp('halfprecision( ) has not been compiled to be multi-threaded.');
   disp('If you have an OpenMP compliant compiler, you should recompile');
   disp('halfprecision.c manually with the appropriate OpenMP compiler');
   disp('flags so that halfprecision( ) can multi-thread. See the file');
   disp('halfprecision.m for instructions on how to do this.');
   disp(' ');
end
end

%--------------------------------------------------------------------------

function halfprecision_test_2019a
disp('-----------------------------------------------------------------------');
disp(' ');
disp('SPOT CHECKS ...');

disp(' ')
disp('Single Value Tests')
d=[1059532000   1059996639   1020952788   1054452098   1019540317   1051411838   1057717945,...
   1057321013   1035188900   1062871840   1036453291   1048400589   1047333619   1041169316,...
   1062279124   1055725692   1065338969   1059140187   1053016952   1052966526   1049355082,...
   1047466717   1059960811   1057962687   1062550723   1043070116   1055788239   1057958540,...
   1054759542   1055267381   1055839636   1062881681   1064978966   1029651805   1059590217,...
   1043139985   1063527922   1061083791   1058967701   1043499464   1061178900   1032632487,...
   1054923658   1060697125   1062044831   1059727808   1032588739   1042488201   1061746082,...
   1026839353   1034904536   1034326491   1054296646   1061179672   1063817555   1048668115,...
   1038345421   1050080045   1025898176   1044597996   1060401766   1060315199   1032671948,...
   1064956484   1054968045   1050786689   1062330574   1051170067   1050004508   1053860464,...
   1050575638   1052733729   1040441505   1062313436   1058902427   1055762301   1061291277,...
   1051381267   1062582256   1061784222   1048567108   1061845687   1057090429   1055478227,...
   1042010389   1037289238   1058349339   1059387663   1064011193   1033344482   1057816651,...
   1055620734   1061517167   1064488589   1033338452   1062238190   1040303873   1061549039,...
   1064020757   1024390200   1053150692   1064960364   1064590360   1046584671   1058715164];
s = typecast(uint32(d),'single');
f = [single([inf -inf nan -nan 0]) s];
hp = halfprecision(f);
hu=[31744   64512   65024   32256       0   14649   14706    9940   14029    9768   13658   14428,...
    14380   11678   15057   11832   13291   13160   12408   14985   14185   15358   14602   13854,...
    13848   13407   13177   14702   14458   15018   12640   14192   14457   14067   14129   14199,...
    15058   15314   11002   14657   12648   15137   14839   14581   12692   14850   11366   14087,...
    14792   14956   14673   11360   12569   14920   10659   11643   11573   14010   14851   15173,...
    13323   12063   13496   10544   12826   14756   14745   11371   15312   14092   13582   14991,...
    13629   13486   13957   13556   13820   12319   14989   14573   14189   14864   13654   15022,...
    14924   13311   14932   14351   14155   12511   11934   14505   14632   15196   11453   14440,...
    14172   14892   15254   11452   14980   12302   14896   15197   10360   13870   15312   15267,...
    13069   14550];
if( isequal(hp,hu) )
    disp('Single Spot Checks Passed.');
else
    disp('Single Spot Checks DID NOT PASS!  Contact author ...');
end

disp(' ')
disp('Double Value Tests')
d=[3568336280   1068636534   1182794464   1066821489   3346998770   1070688574    972352152,...
   1068692904   2511240440   1072493290   2186047606   1072480567    512087016   1070790402,...
   3223546758   1071666481   2422002862   1070892392   2201376514   1071479588   2199615482,...
   1071788983   1124076039   1072374209    806939196   1069923749   1048173454   1071187733,...
     52308343   1072186756     94838997   1072064833   2091863619   1072561441   1798696092,...
   1071507719     44405498   1071579714   4232754981   1072256896    940961487   1072167010,...
   3431564792   1072477558     22540336   1068784596   1499213671   1072343884   1104009292,...
   1071746858   2608866717   1071996448   3203919378   1071452543   1011532312   1072422571,...
   3461353166   1072397335   2810259285   1072634114   1432420856   1071843097   3483663235,...
   1071790717    691025374   1071356641   3824997100   1070446580    190301986   1070987656,...
   3945904972   1072338930   3300823200   1072634577   1314067136   1071995256   1240808464,...
   1069991158    662052244   1072335540   1108664596   1069602027    393125388   1070245655,...
   1138527864   1070686459    502216952   1071577477   3468163488   1072369541   2076190986,...
   1070769715   2797587826   1071216825   2278986706   1071246836    905106952   1069684641,...
     46076940   1070728788     23979171   1072436158   3697698728   1071351943   3462956064,...
   1067474662   3837110105   1072471020   2868318606   1071072066   2203337216   1067866463,...
   4017188175   1072627331   1633575648   1070010557   1787503974   1071276926   1502386650,...
   1071127226   2013206208   1071942680    129069826   1071593474   3155371916   1069897833,...
   2426204506   1071534127   2987824434   1072345631   1372702750   1070715065   1581059218,...
   1070728885   2704735344   1068386823    105348310   1070711469   3945119878   1071380473,...
   3033028554   1071018278   2488982468   1072213944    127866865   1071869929   1731582209,...
   1072518596   1761732568   1068864188   2285776320   1070739766    522194074   1071946804,...
   3828373196   1071942412    999983034   1070951741    312022165   1072636170   1125985304,...
   1070345444   2919087159   1072675919   3555866298   1072601357   4200378168   1070277538,...
   1034399897   1071936439   3409105911   1072523600    135124692   1069642903   2749930620,...
   1071430486   3445242260   1072424553   2896530495   1072628955   2553478916   1070393474,...
   2549743872   1064915826   1111287198   1072658098   4183648000   1071500377   1541436588,...
   1072674252   2444029045   1072182408   1597249920   1068398375   4170022854   1071327395,...
    721453936   1072110287   1074358472   1071214047];
s = typecast(uint32(d),'double');
f = [inf -inf nan -nan 0 s];
hp = halfprecision(f);
hu=[31744   64512   65024   32256       0   11398    9626   13402   11453   15165   15152   13502,...
    14357   13601   14175   14477   15048   12655   13890   14865   14746   15231   14202   14273,...
    14934   14846   15149   11543   15019   14436   14680   14148   15096   15071   15302   14530,...
    14479   14055   13166   13694   15014   15303   14678   12721   15011   12341   12970   13400,...
    14270   15044   13482   13918   13947   12422   13442   15109   14050   10264   15143   13777,...
    10646   15296   12740   13977   13831   14627   14286   12630   14228   15021   13428   13442,...
    11155   13425   14078   13724   14892   14556   15189   11621   13452   14631   14627   13659,...
    15304   13067   15343   15270   13001   14621   15194   12381   14127   15098   15297   13114,...
     7765   15326   14195   15341   14861   11166   14026   14791   13915];
if( isequal(hp,hu) )
    disp('Double Spot Checks Passed.');
else
    disp('Double Spot Checks DID NOT PASS!  Contact author ...');
end

end

%--------------------------------------------------------------------------

function halfprecision_test_2019b
disp('-----------------------------------------------------------------------');
disp(' ');
disp('FULL TESTS ...');
% Load the half( ) function into memory
half(0);
n = 0;

disp(' ')
disp('Large Range Single Tests')
f = randn(1,10000000,'single')*2^(15);
n = n + halfprecision_compare(f);

disp(' ')
disp('Small Denormal Single Tests')
f = randn(1,10000000,'single')*2^(-14);
n = n + halfprecision_compare(f);

disp(' ')
disp('-NaN/Inf Single Tests')
f = typecast(uint32(4286578688):uint32(4294967295),'single');
n = n + halfprecision_compare(f);

disp(' ')
disp('+NaN/Inf Single Tests')
f = typecast(uint32(2139095040):uint32(2147483647),'single');
n = n + halfprecision_compare(f);

disp(' ');
disp('Random bit pattern Single Complex Tests')
umax = double(uint32(realmax));
f = complex(typecast(uint32(ceil(rand(1,10000000)*umax)),'single'), ...
            typecast(uint32(ceil(rand(1,10000000)*umax)),'single'));
n = n + halfprecision_compare(f);

disp(' ')
disp('Large Range Double Tests')
f = randn(1,10000000)*2^(15);
n = n + halfprecision_compare(f);

disp(' ')
disp('Small Denormal Double Tests')
f = randn(1,10000000)*2^(-14);
n = n + halfprecision_compare(f);

disp(' ');
disp('Random bit pattern Double Complex Tests')
umax = double(uint32(realmax));
f = complex(typecast(uint32(ceil(rand(1,20000000)*umax)),'double'), ...
            typecast(uint32(ceil(rand(1,20000000)*umax)),'double'));
n = n + halfprecision_compare(f);

disp(' ');
disp('Random bit pattern eps Tests')
p = -15:15;
t = 2.^p;
d = randn(1,numel(t)+5).*[inf -inf nan -nan 0 t];
r = randn(1,10000000) .* 2.^(round(rand(1,10000000)*32-16));
f = [d r];
hm = half(f);
n = n + halfprecision_compare_fun(hm,@eps,'eps');

disp(' ');
disp('Random bit pattern isinf Tests')
f = randn(1,10000000,'single')*2^(15);
z = rand(size(f)) < 0.3;
f = f ./ z;
hm = half(f);
n = n + halfprecision_compare_fun(hm,@isinf,'isinf');

disp(' ');
disp('Random bit pattern isnan Tests')
f = rand(1,10000000);
f(f<0.6) = 0;
z = rand(size(f)) < 0.6;
f = f ./ z;
hm = half(f);
n = n + halfprecision_compare_fun(hm,@isnan,'isnan');

disp(' ')
fprintf('Total number of mismatches = %d\n\n',n);
end

%--------------------------------------------------------------------------

function n = halfprecision_compare(s)
n = 0;
c = class(s);
fprintf('MATLAB %s -> half timing: ',c)
tic
hm = half(s);
toc
hi = storedInteger(hm);
fprintf('Mex %s -> half timing:    ',c)
tic
hp = halfprecision(s);
toc
fprintf('MATLAB half -> %s timing: ',c)
tic
sm = cast(hm,c);
toc
fprintf('Mex half -> %s timing:    ',c)
tic
sp = halfprecision(hp,c);
toc
if( isequal(hi,hp) )
    fprintf('%s --> Half Precision values match\n',c)
else
    n = n + sum(hi~=hp);
    fprintf('%s --> Half Precision values DO NOT MATCH!!!  (%d of them)\n',c,sum(hi~=hp))
end
if( isreal(sm) )
    usm = typecast(sm,'uint32');
    usp = typecast(sp,'uint32');
    ie = isequal(usm,usp);
    d = sum(usm~=usp);
else
    usmr = typecast(real(sm),'uint32');
    uspr = typecast(real(sp),'uint32');
    usmi = typecast(imag(sm),'uint32');
    uspi = typecast(imag(sp),'uint32');
    ie = isequal(usmr,uspr) && isequal(usmi,uspi);
    d = sum(usmr~=uspr) + sum(usmi~=uspi);
end
if( ie )
    fprintf('Half Precision --> %s values match\n',c)
else
    n = n + sum(usm~=usp);
    fprintf('Half Precision --> %s values DO NOT MATCH!!! (%d of them)\n',c,d)
end
ninf = sum(halfprecision(hp,'isinf'));
nnan = sum(halfprecision(hp,'isnan'));
nnormal = sum(halfprecision(hp,'isnormal'));
ndenormal = sum(halfprecision(hp,'isdenormal'));
fprintf('Total number of values tested = %d\n',numel(s));
fprintf('Number of each kind of value tested:\n');
fprintf('inf = %d , nan = %d , normal = %d , denormal = %d\n',ninf,nnan,nnormal,ndenormal);
end

%--------------------------------------------------------------------------

function n = halfprecision_compare_fun(hm,fun,funstr)
n = 0;
hi = storedInteger(hm); hp = hi;
fprintf('MATLAB %s() timing: ',funstr)
tic
hmf = fun(hm);
toc
fprintf('Mex %s() timing:    ',funstr)
tic
hpf = halfprecision(hi,funstr);
toc
if( isequal(class(hpf),'uint16') )
    hmf = storedInteger(hmf);
end
if( isequal(hmf,hpf) )
    fprintf('Half Precision %s() values match\n',funstr)
else
    n = n + sum(hmf~=hpf);
    fprintf('Half Precision %s() values DO NOT MATCH!!!  (%d of them)\n',funstr,sum(hmf~=hpf))
end
ninf = sum(halfprecision(hp,'isinf'));
nnan = sum(halfprecision(hp,'isnan'));
nnormal = sum(halfprecision(hp,'isnormal'));
ndenormal = sum(halfprecision(hp,'isdenormal'));
fprintf('Total number of values tested = %d\n',numel(hm));
fprintf('Number of each kind of value tested:\n');
fprintf('inf = %d , nan = %d , normal = %d , denormal = %d\n',ninf,nnan,nnormal,ndenormal);
end
