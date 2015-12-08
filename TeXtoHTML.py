#! /usr/bin/python

"""
  ooo        ooooo           oooooooooo.             ooooo      ooo
  `88.       .888'           `888'   `Y8b            `888b.     `8'
   888b     d'888   .ooooo.   888      888  .ooooo.   88 `88b.   8   .oooo.
   8 Y88. .P  888  d88' `88b  888      888 d88' `88b  88   `88b. 8  `P  )88b
   8  `888'   888  888   888  888      888 888ooo888  88     `88b8   .oP"888
   8    Y     888  888   888  888     d88' 888    .o  88       `88  d8(  888
  o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o88o        8  `Y888""8o

@summary:      Script for converting a LaTeX document to various utput formats.
               The script does **not** work on a Windows computer.


@author:       Sigve Karolius
@organization: Department of Chemical Engineering, NTNU, Norway
@contact:      sigveka@ntnu.no
@license:      GPLv3
@requires:     * Python 2.7.6 (python --version)
               * Pandoc 1.12.4.2 (pandoc --version)
                 >* pandoc-citeproc
                 >* pandoc-crossref
               * ImageMagick 6.9.0-0 (convert --version)
@since:        15.06.2015 (SK)
@version:      0.1
@todo 1.0:
@change:       started
@note:         Work in progress.

port install ghc, hs-cabal-install
cabal update
cabal install pandoc, pandoc-citeproc, pandoc-crossref
echo "export PATH=${PATH}:~/.cabal/bin/" >> ~/.profile


apt-get install ghc, cabal-install
cabal update
cabal install pandoc, pandoc-citeproc, pandoc-crossref
echo "export PATH=${PATH}:~/.cabal/bin/" >> ~/.bashrc

"""

import os, sys, re



files = []

with open('Main.tex','r') as FILE:
  for line in FILE:
    match= re.match(r'\\input\{(.*)\}|\\include\{(.*)\}',line)
    if match:
      files.append(match.group(1))

with open('tmp.tex','w') as tmp:
  for f in files:
    with open('%s' %(f)) as FILE:
      for line in FILE:
        fig=re.match(r'^\s*\\includegraphics.*\{(.*)(\.eps|\.pdf)\}',line)
        if fig:
          os.system("convert %s%s %s.png" %(fig.group(1),fig.group(2),fig.group(1)))
          line = re.sub(r'%s' %(fig.group(2)), '.png', line)

        tmp.write(line)

#os.system("pandoc --standalone --from=latex --toc --filter pandoc-crossref --bibliography=./Content/Bibliography/bibliography.bib --highlight-style tango --to=markdown_github --output=tmp.md tmp.tex")

os.system("pandoc --standalone --from=latex --bibliography=./Content/Bibliography/bibliography.bib --to=markdown_github --output=tmp.md tmp.tex")
os.remove('tmp.tex')


with open('edit.md','w') as newFILE:
  with open('tmp.md','r') as oldFILE:

    for line in oldFILE:

      if re.match(r'^!.*data-label.*', line):
        lbl = re.match(r'.*fig:(.*)".*', line)
        line = re.sub(r'\n' ,'{#fig:%s}\n' %lbl.group(1), line)

      while re.match(r'.*\\\[[A-Za-z]+:[A-Za-z0-9]+\\\].*',line):
        ref = re.match(r'.*(\\\[([A-Za-z]+):([A-Za-z0-9]+)\\\]).*',line)
        line = re.sub(r'%s' %'\\\[[A-Za-z]+:[A-Za-z0-9]+\\\]', '[@%s:%s]' %(ref.group(2),ref.group(3)), line)

      eqStart=re.match(r'.*(\\\[).*',line)
      eqStop =re.match(r'.*(\\\]).*',line)
      
      while re.match(r'^\n{0,1}.*(\\\[).*',line):
        pos=eqStart.start()
        line = re.sub(r'\\\[', r'\n\n\[', line)

      if re.match(r'.*(\\\]).*',line):
        pos = eqStop.start()
        line = re.sub(r'\\\]', r'\]\n\n', line)

#      if eqStart:
#        pos=eqStart.start()
#        line = re.sub(r'\\\[', r'\n\n\[', line)
#
#      if eqStop:
#        pos = eqStop.start()
#        line = re.sub(r'\\\]', r'\]\n\n', line)

      newFILE.write(line)



os.remove("tmp.md")

os.system("pandoc --standalone --number-sections --smart --self-contained --toc --filter pandoc-crossref --highlight-style tango --mathml --from=markdown_github --to=html --output=Main.html edit.md")

os.system("pandoc --standalone --smart --toc --filter pandoc-crossref --bibliography=./Content/Bibliography/bibliography.bib --highlight-style tango --from=markdown_github --to=docx --output=Main.docx edit.md")

#os.remove("edit.md")

