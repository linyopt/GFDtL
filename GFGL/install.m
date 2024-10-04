%% Install script
% Checks dependencies for GFGL ie SLEP and GFLseg are accessible to MATLAB
% path...if not download them
% Alex Gibberd 2015

% This part gets the unix path and adds the local directory...this was
% required to work on mac. It may be irrelevent for your system setup.
PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin']);

download_dir=[pwd,'/libs/'];

% Check for GFLSeg folder/file
if (~(exist('GFLseg-1.0','dir') && exist('gflasso','file')))  
    display('GFLSeg not found on path...attempting to download');
    % GFLseg package from Bleakley & Vert's webpage
    url='http://members.cbio.mines-paristech.fr/~jvert/svn/GFLseg/html/GFLseg-1.0.tar.gz';
    target=[download_dir,'GFLseg-1.0.tar.gz'];
    % IF OSX use
    unix(['curl -o ' target ' ' url ]);
    
    % If unix OS use below
%     unix(['wget -O ' target ' --timeout=100 "' url '"']);
    %[s,r]=unix(['tar -C  -zxvf GFLseg-1.0.tar.gz']);
    display(pwd);
    [s,r]=unix(['tar -zxvf libs/GFLseg-1.0.tar.gz -C libs/']);
    if(s==0)
        display('installed GFLSeg successfully');
    else
        display('There was an error installing GFLSeg');
        return;
    end   
end

if(~exist('fusedLeastR','file'))
    % Download SLEP
    display('Required slep files not found...attempting to download');
    % From Sparse Learning with Efficient Projections website
    url='http://yelab.net/software/SLEP/SLEP_package_4.1.zip';
    target=[download_dir,'SLEP_package_4.1.zip'];
    unix(['curl -o ' target ' ' url ]);

%     unix(['wget -O ' target ' --timeout=100 "' url '"']);
    display(pwd)
    [s,r]=unix(['unzip libs/SLEP_package_4.1.zip -d libs/']);
    if(s==0)
        display('installed SLEP successfully');
    else
        display('There was an error installing SLEP');
        return;
    end   
    display('MEXING SLEP');
    run('libs/SLEP_package_4.1/mexC.m'); 
    
end

% update path
addpath(genpath(pwd));

% Check dependencies again
if(~(exist('fusedLeastR','file') && exist('GFLseg-1.0','dir') && exist('gflasso','file')))
   display('Instalation failed..please try downloading by hand'); 
else
   display('Instalation succesful'); 
end
