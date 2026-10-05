#!/usr/bin/env python3
"""Build the pinned Gmsh/OCCT runtime inside this repository.

Offline by default. Put the two verified archives in build/native-deps/downloads,
or pass --archives-dir. --download explicitly permits fetching missing archives.
On macOS the selected compiler is Apple's /usr/bin/clang++; Linux defaults to c++.
The Linux configuration is supplied for portability, not a Linux verification claim.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / 'dependencies.lock.json'
OCC_LIBRARIES = ('TKDESTEP TKDEIGES TKXSBase TKOffset TKFeat TKFillet TKBool TKMesh '
                 'TKHLR TKBO TKPrim TKShHealing TKTopAlgo TKGeomAlgo TKBRep TKGeomBase '
                 'TKG3d TKG2d TKMath TKernel').split()

def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()

def checked_archive(item, archives, downloads, allow_download):
    target = downloads / item['archive']
    source = next((base/item['archive'] for base in [downloads, archives]
                   if base and (base/item['archive']).is_file()), None)
    if source is None:
        if not allow_download:
            raise RuntimeError('Missing '+item['archive']+'; use --archives-dir or explicitly --download')
        downloads.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(target.suffix+'.partial')
        with urllib.request.urlopen(item['url']) as response, temporary.open('wb') as stream:
            shutil.copyfileobj(response, stream)
        if digest(temporary) != item['sha256']:
            raise RuntimeError('Downloaded archive checksum mismatch: '+str(temporary))
        temporary.replace(target)
        source = target
    if digest(source) != item['sha256']:
        raise RuntimeError('Archive checksum mismatch: '+str(source))
    return source

def extract_or_verify(item, archive, sources):
    destination = sources/item['source_directory']
    with tarfile.open(archive) as tf:
        members = tf.getmembers()
        for member in members:
            target = (sources/member.name).resolve()
            if not target.is_relative_to(sources.resolve()):
                raise RuntimeError('Unsafe archive entry: '+member.name)
            if member.issym() or member.islnk():
                link = ((target.parent if member.issym() else sources)/member.linkname).resolve()
                if not link.is_relative_to(sources.resolve()):
                    raise RuntimeError('Unsafe archive link: '+member.name)
        if destination.exists():
            ignored = {item['source_directory']+'/'+s for s in item.get('generated_cache_paths',[])}
            for member in members:
                if member.isfile() and member.name not in ignored:
                    target=sources/member.name
                    if not target.is_file() or digest(target)!=hashlib.sha256(tf.extractfile(member).read()).hexdigest():
                        raise RuntimeError('Existing dependency source differs from archive: '+str(target))
        else:
            sources.mkdir(parents=True,exist_ok=True)
            if sys.version_info >= (3,12):
                tf.extractall(sources,filter='data')
            else:
                tf.extractall(sources)

def build_commands(args, root):
    prefix=root/'install'
    mac=platform.system()=='Darwin'
    if platform.system() not in ('Darwin','Linux'):
        raise RuntimeError('This bootstrap supports macOS and Linux only')
    cc=args.cc or ('/usr/bin/clang' if mac else shutil.which('cc'))
    cxx=args.cxx or ('/usr/bin/clang++' if mac else shutil.which('c++'))
    if not cc or not cxx: raise RuntimeError('C and C++ compilers are required')
    cmake=shutil.which('cmake')
    if not cmake: raise RuntimeError('cmake was not found')
    rpath=('@loader_path' if mac else '$ORIGIN')+';'+str(prefix/'lib')
    common=['-G','Unix Makefiles','-DCMAKE_BUILD_TYPE=Release','-DCMAKE_C_COMPILER='+cc,
            '-DCMAKE_CXX_COMPILER='+cxx,'-DCMAKE_CXX_STANDARD=17',
            '-DCMAKE_INSTALL_PREFIX='+str(prefix),'-DCMAKE_INSTALL_RPATH='+rpath,
            '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF']
    if mac:
        common+=['-DCMAKE_OSX_ARCHITECTURES='+platform.machine(),
                 '-DCMAKE_OSX_DEPLOYMENT_TARGET='+args.macos_deployment_target,'-DCMAKE_MACOSX_RPATH=ON']
    ob=root/'occt-build';gb=root/'gmsh-build'
    occ=['-DBUILD_LIBRARY_TYPE=Shared','-DBUILD_CPP_STANDARD=C++17','-DBUILD_RELEASE_DISABLE_EXCEPTIONS=OFF',
         '-DBUILD_MODULE_FoundationClasses=ON','-DBUILD_MODULE_ModelingData=ON','-DBUILD_MODULE_ModelingAlgorithms=ON',
         '-DBUILD_MODULE_ApplicationFramework=ON','-DBUILD_MODULE_Visualization=OFF','-DBUILD_MODULE_DataExchange=OFF',
         '-DBUILD_ADDITIONAL_TOOLKITS=TKDESTEP TKDEIGES TKXCAF TKVCAF',
         '-DBUILD_MODULE_Draw=OFF','-DBUILD_MODULE_DETools=OFF','-DBUILD_Inspector=OFF','-DBUILD_SAMPLES_QT=OFF',
         '-DBUILD_USE_PCH=OFF','-DBUILD_RESOURCES=OFF','-DINSTALL_TEST_CASES=OFF','-DINSTALL_DIR='+str(prefix),
         '-DUSE_TBB=OFF','-DUSE_FREETYPE=OFF','-DUSE_FREEIMAGE=OFF','-DUSE_FFMPEG=OFF','-DUSE_OPENVR=OFF','-DUSE_XLIB=OFF']
    suffix='.dylib' if mac else '.so'
    libs=';'.join(str(prefix/'lib'/('lib'+name+suffix)) for name in OCC_LIBRARIES)
    search=[str(prefix)]
    freetype=args.freetype_prefix
    if freetype is None and mac and Path('/opt/homebrew/opt/freetype').is_dir():
        freetype=Path('/opt/homebrew/opt/freetype')
    if freetype: search.append(str(freetype.resolve()))
    gmsh=['-DGMSH_RELEASE=ON','-DENABLE_BUILD_DYNAMIC=ON','-DENABLE_BUILD_LIB=OFF','-DENABLE_BUILD_SHARED=ON',
          '-DENABLE_OCC=ON','-DENABLE_OCC_CAF=ON','-DENABLE_OCC_STATIC=OFF',
          '-DOCC_INC='+str(prefix/'include/opencascade'),'-DOCC_LIBS='+libs,
          '-DCMAKE_PREFIX_PATH='+';'.join(search),'-DCMAKE_LIBRARY_PATH='+str(prefix/'lib'),
          '-DCMAKE_INCLUDE_PATH='+str(prefix/'include'),'-DENABLE_FLTK=OFF','-DENABLE_GRAPHICS=OFF',
          '-DENABLE_MPI=OFF','-DENABLE_OPENMP=OFF','-DENABLE_OCC_TBB=OFF','-DENABLE_WRAP_PYTHON=OFF',
          '-DENABLE_WRAP_JAVA=OFF','-DENABLE_TESTS=OFF','-DENABLE_PETSC=OFF','-DENABLE_MUMPS=OFF',
          '-DENABLE_CGNS=OFF','-DENABLE_MED=OFF','-DENABLE_CAIRO=OFF','-DENABLE_GMP=OFF',
          '-DENABLE_NETGEN=OFF','-DENABLE_SYSTEM_CONTRIB=OFF','-DENABLE_P4EST=OFF','-DENABLE_POPPLER=OFF']
    return [('occt-configure',[cmake,'-S',str(root/'src/OCCT-7_9_3'),'-B',str(ob)]+common+occ),
            ('occt-build',[cmake,'--build',str(ob),'--parallel',str(args.jobs)]),
            ('occt-install',[cmake,'--install',str(ob)]),
            ('gmsh-configure',[cmake,'-S',str(root/'src/gmsh-4.12.2-source'),'-B',str(gb)]+common+gmsh),
            ('gmsh-build',[cmake,'--build',str(gb),'--parallel',str(args.jobs)]),
            ('gmsh-install',[cmake,'--install',str(gb)])]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-root',type=Path,default=REPO.parent/'build/deps/cad')
    parser.add_argument('--archives-dir',type=Path)
    parser.add_argument('--download',action='store_true')
    parser.add_argument('--jobs',type=int,default=8)
    parser.add_argument('--cc');parser.add_argument('--cxx')
    parser.add_argument('--freetype-prefix',type=Path)
    parser.add_argument('--macos-deployment-target',default='14.7')
    parser.add_argument('--print-plan',action='store_true',help='Print commands without modifying any files')
    parser.add_argument('--verify-archives-only',action='store_true')
    args=parser.parse_args()
    if not 1<=args.jobs<=256: parser.error('--jobs must be in 1..256')
    root=args.build_root.resolve();commands=build_commands(args,root)
    if args.print_plan:
        print(json.dumps(commands,indent=2));return
    manifest=json.loads(MANIFEST.read_text())
    archives={item['name']:checked_archive(item,args.archives_dir,root/'downloads',args.download)
              for item in manifest['dependencies']}
    if args.verify_archives_only:
        print('Verified both pinned dependency archives.');return
    root.mkdir(parents=True,exist_ok=True);logs=root/'logs';logs.mkdir(exist_ok=True)
    for item in manifest['dependencies']:
        archive=archives[item['name']];target=root/'downloads'/item['archive']
        target.parent.mkdir(exist_ok=True)
        if archive.resolve()!=target.resolve():shutil.copy2(archive,target)
        extract_or_verify(item,target,root/'src')
    state={'schema':'v3-native-dependency-build-1','status':'running','manifest_sha256':digest(MANIFEST),
           'prefix':str(root/'install'),'platform':platform.platform(),'commands':[],
           'sources':manifest['dependencies'],'linux_verification':'unperformed'}
    report=logs/'provenance.json'
    env=dict(os.environ,LC_ALL='C',LC_CTYPE='C',LANG='C')
    for key in ['CASROOT','OpenCASCADE_DIR','OCC_DIR','DYLD_LIBRARY_PATH','DYLD_FALLBACK_LIBRARY_PATH']:
        env.pop(key,None)
    try:
        for label,argv in commands:
            print(label,flush=True);entry={'label':label,'argv':argv,'started':time.time()};state['commands'].append(entry)
            report.write_text(json.dumps(state,indent=2)+'\n')
            with (logs/(label+'.log')).open('w') as stream:
                proc=subprocess.run(argv,env=env,stdout=stream,stderr=subprocess.STDOUT)
            entry.update(return_code=proc.returncode,elapsed_seconds=time.time()-entry['started'])
            if proc.returncode:raise RuntimeError(label+' failed; see '+str(logs/(label+'.log')))
            if label=='gmsh-configure':
                text=(root/'gmsh-build/src/common/GmshConfig.h').read_text()
                if '#define HAVE_OCC' not in text or '#define HAVE_OCC_CAF' not in text:
                    raise RuntimeError('Gmsh is missing OCC/CAF; verify installed Freetype development files')
        version=subprocess.run([str(root/'install/bin/gmsh'),'-version'],env=env,text=True,capture_output=True,check=True)
        if (version.stdout+version.stderr).strip()!='4.12.2':raise RuntimeError('Unexpected Gmsh version')
        state['status']='built_requires_runtime_verification'
        print('Built '+str(root/'install')+'. Run the project native/MATLAB runtime tests before installation.')
    except Exception as error:
        state['status']='failed';state['error']=str(error);raise
    finally:report.write_text(json.dumps(state,indent=2)+'\n')

if __name__=='__main__':main()
