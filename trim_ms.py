t=open('manuscript.tex').read()
pos=t.find(r'\section*{Supplementary Information}')
doc_end=t.find(r'\end{document}')
t=t[:pos-1]+'\n'+t[doc_end:]
open('manuscript.tex','w').write(t)
print(f'Trimmed to {len(t)} chars')
