from pathlib import Path
from PIL import Image, ImageDraw
folder=Path(r'E:\脑膜瘤npj\tmp\manuscript_render')
files=sorted(folder.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[-1]))
thumbw=220; thumbs=[]
for i,f in enumerate(files,1):
    im=Image.open(f).convert('RGB')
    h=round(im.height*thumbw/im.width)
    t=im.resize((thumbw,h))
    card=Image.new('RGB',(thumbw,h+22),'white'); card.paste(t,(0,22))
    ImageDraw.Draw(card).text((5,4),f'Page {i}',fill='black')
    thumbs.append(card)
cols=4; rows=(len(thumbs)+cols-1)//cols; cellh=max(im.height for im in thumbs)
sheet=Image.new('RGB',(cols*thumbw,rows*cellh),(225,225,225))
for i,im in enumerate(thumbs): sheet.paste(im,((i%cols)*thumbw,(i//cols)*cellh))
sheet.save(folder/'contact.png')
print(f'{len(files)} pages; contact sheet saved')
