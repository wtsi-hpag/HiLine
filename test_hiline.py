# Copyright (c) 2020 Ed Harry, Wellcome Sanger Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import bz2
import os
import re
from subprocess import Popen, PIPE, STDOUT
from threading import Thread

from HiLine import Pipeline
from HiLine import version

reads = b"BZh91AY&SY9\x9c\xd8)\x00\x01\xd1\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xe0\x06\xbf\x03\xd7\xdd\xd1\xea\xad{|\xe7\xbe\x8e\xf7\xb9\xdd\xe6\xe7{\xcb\xb6\xf8\xca\x9ePh\xf2\x8d\xa9\x91\xa3O$\xc4\xc0&5<S2CM\x07\xa0F\t\x8dF\x9e\x89\xe9\x1bS\x13M\x01\x88\xf4\xd4f\xa6\x08\xc8\xcd'\xa6I\xa7\xa9\xa1\x90\x1a4\xd1\xa7\xa4\xf2hj4\xd3 f\x84\xc3S\xca=2jmG\xa0\xd46I\x93\xd4\xf5\x06G\xa8eL\x04\xc4h\xf4M<\xa1\xeazO(\xf4M\x1e\x91\xea4f\xa0\xc6\xa0\xc9\xa3&\xc9\xa2z\x9eS'\xa4\xd3F\x9927\xa94\xf4\x135\x1e\xa6\x994\xf2\x8c\x9a4\xc4l\x8d4Fh\x9bD\xf2\x9ag\xa4\x9b$\xf4\xc9\r\xa9\xe9\x03L\xca\r\x06CC\x11\xa7\xa4\xc2\x1a\x04\x8c\x8c'\xa9\x80jd\xd3\x1a\x1a\t\x88\xc2hh\xc4\xc0\x02z\x04a\xa1\x06\x1a54\xfdPd`\x18A<\x82x\x9a4&\x06\xa350\x00\x13OH\xcd'\xa8\xd3#LF\x98!\x80\xd0\t\x80\xd10#\x02\x13\xda\xa7\x89=Fd\x9e\x93F&\x13jmM\x8ah3H\xd3\xd0h\xd4zh\xc8\x9a=F\x9bG\xa54i\xeaoQ<\x9a\x9abdcB\x0fP\xf5\x18L\x9e\xa6\x11\xea\x1a=C@\xf50\x83\t\xb4\x80\x1a0\x98G\xa9\x93C\xd04\x8fQ\xa6#i0\x8c\x8c:\x8fQ\xea=G\xa9\xea=#zQ\xa3\xd3S\xc53\x11\xe2\xa6\x9e\x91\xfa\xa7\xeaCi\r4\x1e\xa3\xf55\x1e\xa6\x83\x11\xa1\xb5\r2z\x87\xa4\xd0\x07\xa6\x9a\x9b(\xf5=Ci\x0fS\xd444h\x00\x0f@\xd2a\x03\xd47\xaa\x19\x1a\r\xa8i\xeaz\x13M\x19?T\xf5\x03\xd4\xf5\x06\x8d2\x04\x14\x10MO\t=&j2a\xa4\xc4\xf5\x01\x88\xf52i\xa6M\x0fF\xa6F\x800@\x0c\x9e\xa0\x0f(\xf54\x0fH\x1e\xa0\xd1\xa0=#\xd4\x06\x9a\x00h\x19\x00\x1e\xa6\x86\x80\xd1\xa0\xc8\xd04\r\x00\x00\x03@4\x8b$b\x06\xe2\\\xda\x1e \x0eE\xe4\xb8\x1eh\x80=\xa22\xb0jS@\xe2%\x06\x18Wc\x04'\x08\xd1\x04ALeI\xd0\x84j\x92\x02\n\xb1\x13\xfe\x8fG\x9e\xb8\xc4 D$\x048\x12\x94\xc3\x18\n`\xd0E\x11\xc4\x8c\xc4\xc6\xcd?~\xcf8fr\xb9\xb9*\xf3\xc0\x82\x01\x9c*\xe0&1\xcc\x8d!\x83,\x14\x04\x85\xbd\"D\x02\x05\xf5\x7f\xa1\x95t\x10@\x1f\xa9Z\x06v\x16\x16H&\xb7\xe7X\x1e\xf4-\x03\x04B\x8cT<\xd6]\xaa\x8c\x10\x82l\xc7\x8e\xde\xd1)\xbf\xd5\xdfR%d\x95\x93\xbb\xe2\xc3\x16\x01\xd4LOd\xe3\xc5\xd8\x7f\xbaJ\xa8\x0fA\x0bu\xa2/k\x99\xd6\x95\x92\xd2\x8cY\xd6\xbc\xcbG\xe4!\xb6\xfdV\x1e\xfed{\xecQ\x8a\x17\xba\xbbD\x86R\xebf\x84r\xa3\x88\xe2\xc3q\x1a\xa5qR\x8d\\A\xe0\x8b\xf3\xd2\xe0S\xc3\xbf!\xdd\xbc*;\x93:\xef\xb6\x98\x9f\x10\xbf\x9c3Q\xca)\xb5\xf2\\\xaaY\xca\xefq\xa2v\xd4\xea\x85\x91\xfeG\xec\xf9\x8a\xe3D\x87\xdb@\x1f\xcf\xa1\xc2\xbf\xcdc\x90\xb6/\xdd\x9c&\xb6d]h/\x88\xd7j\x11\x8c\xe9N\x9f4\x1cg\x8b\xa1\xa4\x8c$%v\xb5T\x8d:\xf4\x14\x913s))!{8\xea\xb7-\xd1\x08:\xa7\xf9\xaa*V\xed'\x91}\xb2\xed\xf89\x19g\xa9\x10\xf1\xfaS\xe8yxT3\xd01,\x8a\x8a\xed\xe0P\xfe=\x88\xd1\x15\x0e[W\xa1\xb2\xc0\xb3\xf8\x8e\x01\xd3\xcaA}\xd6yO8\x95f\x81~\x84\xea\xb9\xc9e\x11JK\xa6\xdd\t\xd6\xa3+\xf7\xf4~W\x19g\x1e\xdf\xc3#\x95\x86\x14\x8dA\x0b.\xc9\x11\x80\xf5\xfc\xd5:\xa0a\x02 \xc4\xe7<\xcc&\x8f\xfd\xd9\x887H\x03\xc1\xae\xfb\xf2\xdb\xd7\xf8\x1e\xb4\xb2\xf2\x8c\xc0\x0e$\x96;9+V\xdb\x0ej\xa0\xfa\x87}\xc4\xd4\xbf\xac\x97\xd5#\xee\xec\x07\x94\xf1@O\x16D\x18x+q\xd5X\xeb\xf6\xeax\xc1(t\xee\xcb\x06\xd8\xfeS\xb2\x06\x00\x9c~%\xcd\xb2\xa8\xa9\xe3\xe3h\xedd\xc4\x19\x00\x14 5\xfc\x93mO>(\x1c]\xbe\x90\xea\x9f#\x18J\x15\xd6joc\x972\xb1\x9f\x13U\xcf\x94\x03]A\x7f\x80\x1a\xd8\x13EzWC\x9ez.,\x1a\xae\xdetr\xda\x91^\xb9\xe6\x91WoS\xa0\xe1j\x13\xb0tF\x9aM\xb5\x00\xdb!\x17\x01\xa6e\x1c\xc1\xf5\x99\x8c0\x9f\xb7\xb0\xb4\x92\xdf\xe5\xff0\xd9\xc5\xc4\xc3\xe5F4\x00\xeb1\x8a\xd9'\xd0%\xe9Z;w\xa7,\x93\xd5$\t\x1cXU6'U\xdc\"\xb3\x89\x9f\xafj\xe7\x03?\xbb\xc3\x8d\xd5A~\xc2\x1b\xfc\xbe\xb4E\x87\x0e\r\xb8\xc7\x9f\xf4n$x/\xbe]\xd7\xd9\xe4}v9\x1d}\xb4\xcc\xd8\x8f\x9aA\xb7\xfc\x033\xffl\xf5\xab\xc7\xd3\xc95\x04\xc9\x0f\x8b7w\xf3\xd4&\x1e\\\xc6Y\xf8\xda\x07\x1d<C\xc1r\x16\x90\x05Z\x03\x0e\xd5\xf9\x9c\x00JB2\x00\x8d\xaa\x18\xa3}\x18%\x04\x9b8\xd9\xc3\xa9j\xcc\xddS\xc4\xc8\xb4 \xb3\x1b\x84A\xc0H=\xe7e)\xcezCd\x16\xb6R\xee\xf0\x7f\x8e\x1d\x9b7\x81Xq\xc7\x97\xf2)\x84\xbepQ\xa8\xa2\xd6\x9c\x86\xbb\xafxXT\xd8\x8e\xfaTLK\n\x90\xbet]\x91\x04\xd4\x15\xb6X\x1e\x8d\xc99k\xc8T-\x85\x82\x1cv\x84\x8f\xa5o\xd0,\xcbkEV\xb2 v\xa9e\x1d\xcb<k\xc7\t\xfaF\xb9\xae?\x18\xec\xe8\x82n\xb37\xc7\x99\xc4\xf4\xed\xa2\xc3\x07\xb9z#\xb6\x17\xa6E`\x9d|<\x0c\xfb\x97\xc6\xd3G\xfbkf=\t\xca\x98\xce\x99d\x9a\x16E\x07M\x1aR\xee\xeby\x02\x9b\x1f\x02\x86\xd8\xa9\x87\"m\xb4\x89}\xa0\xc3\x8d\xf1\xe9g\x177\x01\xdar>i\xeeX0\x08(\x15\x02\x89\xd1\xd1E\xaf\xbar\x9dG\x88U\xbbjWX\xb1\x12.@7\x07%7\x19\xb2\x17/\xc0hI'JF\"\x81\xf7\xafP\x1b\n\x1b\x1ad8D\xccIy\x04\n\xd9\x89\xfd\xca\xb6\x19\xe8\xe8\xb7\xca/;d\xc1\x9d\xa1p\xa9\xafS\x19!7\xc6\x9dH2r\xf9omG\xb3\x82\xd6\x9c\x99H\xe4'/\xa3\r\xedGp\xfcq\xdd\xb4\t\xa2\x8e\xaeG\x05\xe0\x05\"\xffh\xaa\xf4\x10\xe0\x07\x86\x9bv3\x99\x13\x87\x85\xe1W(\xcaL\xd1\xbe\x16\x87\xfe\x14\x9cP\x82K\xba\n\xff\"\xa3!\xd3\x04\n\xa8-\x80\x900\xb3\xbds\xe8 \xead\x1e\t\x17\x06\xa4\x9a3\xfb\xee\x03\xdaR\x86V\xfd\r\xfa\xc7\xb0\x90c\xcfcJ-\xbc\x0b\xe3\xd8zy\x035\xa4H+He\xe3\x92\xb82\\vV0\x8b!Cd\x7f*\xc9\xbe8e7\xb6\xcf,M<\xfb{h\x938\xa5\xd0\xd4{\x1c\xe1\x16EDb\x80\x99%\xd1\xf6r\x90dS1K;\xd0\x02\x86s\xb86{\x1dU\xdd&3\x04\xe6\xce_\xa6\x1dx\xc8\x8c{**\xa2\xbc:\xf2p\xa5s+\xd8\xa8_\x16)OH&\xcaL\xec\x18\x8f3K\xdc\xd9\x00\xfdZ\xd4\x14\xc1CGu\xba\x0cV\x052\x80d>\x0b\x989\x9f\xa3}\xa8\xad\xa3\xe3Qa\xc5\xcc\xba\xf6\xc9\xb5\x11GJ\x9eU\xf1\x1d]Th\x007W\xefS\x07\xaa!\xb1Y\x01\x12\xe3#\xf8\xfe\xa2o\xfbM+k\xd2\x8b\x91\xcf\t\xe5\xbbb\x11\x90\x96\x07#h\xbd~\xd68Cd\x89\x96\xa9\x02Q\xa963\xda\xe1\xc7V|\x8e=\xfcK>!\xaa\x92p\xd8D\x17\xe6\t\x05\xdf\xbd?\xbb\x96X\xf0\xf6\x92\x964\xfc\t\xde\xabt\xc8c\x93\xefr\xd9%X\xb2\xb2\xe1\xf3\x8e\x04FC#\xc3\xeb\xe0!\x02\xe1\x9d\x17 e\xf8v*\xcd\xc6eHP\x99\x93\x1a8n\xc2\x97\xf0L\x96G\xb3G\x9e\xf5GV\xe3\x89\xe4\xa8\xba\xf4\nF\x9d\t#N,}:\xd9\x07r\xf7-\\\xc8\xe8\xef\xaa\"\x91M\x11\x08\xfa\xca\xc3\x97'\xe6_S\xa8/(?\no\xdb\xb0\xc0\x9dA\x7fy\x91\x93\xa6\xe3\xf6\x18h\x9ept)eOf\x00\xb7:^@E\x86\xc9\r=\xd9\x9b\x99\xba\x85\xfb6\xad\xc3#\xd6\x03\x04e\x1a$\x9c\x18\xf0\x86\x86\xf6\xf4h\xee\x7f\x15\xaeno\x83\xd6z\x1c\xa1\xc1:J<h\x00\x07\x0b%Uo2&\xdeX\xc3\xf8P\x10\xc2'ZJ\x84\xddZ*\x0e\xcb\xd8\xa4\xddR\xa3L\x88\x0c\xfa\\2\xe7\x1cko>\xca m\xde\x88\x19\x02\x19\xb9\x89@&\\\r)1\xad\xacI\x0c\x86NB\xb2\xe4\x835\xa1\x81\x9e\xf9\xa1\xf2A)/q*\xfe\x91\x1b\x89u\xf9\xa97N\xff6\x02\xf1YXu\xef\xd8SKk\xb7\xabi\x8a\xd5\x81\xa9\xfa\xd4B\xb1\xe3D\xe0\x9c\xae\\^\xf4O\t\x14\xee\x17\x95\xf0\xd1\x9e1\xb9v\xdch\xb1t\xfc2R\x80\xa3\xd0Z\xca\xd5\x07Uf\xfc\xf8\xc0\x07\x02\xcd\x8f<\xac\xe2\x85\x0c\xa2\xe3\x06\xcfI\xe1P\t'\xe5}X\xfa\xb8\xabM\x02\x9eWF\xf9\x98\xf7\xb2\x92\xe5+\xd5\xe5\\\x11$\xb8\xf3\x19\xb7>\xa1\xd0\xc95>\xa7)\xa1\xbc\x08\x03\xd1\xc5\x89\x96\xce)\xa9M\x9e\xbf\x0b\xf61\x03Hj\xbf\xf6Bz\x01_\x88\xea\x03\xa2\xad\x07\"\xbeX\x13\xd9\xb8)B\xa7gi\xdf\xcdm\xa1ra\x95 \xdeUbN\xe2\xc4\x96\xbb1\x8a\x0erB,\xa8y\xda\x0b\x8f\xa1\xe1\xe9!P=\xe4S\x14\x9b\xa3c\xdeA\xdc\xcb\xc8\xbc\xde\xac\xb2\xf2BQ,\xaa\x00\x81l\x1eI\xa54\xe2\xaa\xd87\x83\xaf\xb1\xe8\x9cWX\xaaa\xe9\x1f\xd8\"\xb9y\xcfE\xda\r\n\xa0\xf6g(9\xd92\x90\x15T&$#\x13o_f\xbd\xf0\xc3v\xf4\x94\xd3\x05\xcd\xdcA\xe1\x81\xed\xdb\x90D{x\x97\x92\xe5,'\xa6\xf4\xc30\xdf5\xa5P\x89\x06\x94\x84k`s\xdb&\xbd\xa2\xf0Re\xa1\xc2\x0e\xab\xbaI\xb2H\x00\xeb\xa0\xf2\x08\xd1\x1f\xbc\x06g\xd2N*\x853\xb1\xe9\xd8iA\x8c\x9c\xd4\xce\x0b\x94D\xa7\x9blPT\xae<\x8c\xd1\xb0m \xda\x8b\xc5\xd5\xe0LQ\xb2\xc77\x13\xe2H3Kt5\xf9\xb3\x0b-\xe9\xff\xc6\xa1\x85\x07QWa8\n\xd9\xeb]\x06R\x8d^\xd2,=\xc3K\xf5\x8eK'\xd1l\xbd\xef\x88g\xfa\xed\xae\xeb+\xa5$<M1;\xc1\xfd\xaaA\xc9\xfdlJ_\xcf\xbbm\x85\x8e\x9c|\xf0\x13\xd4\xaf.m\x85^E\xf5s\x0eZ\x8ahna\x11\x14\x91\xf6\x97N\xe8\n_kZ\x98\xce\x9c\x95\xffB\xd8\xc6\x8d\x8b\xb8\xd7u\xc5O\xb8\x05\x9e\xac\x10\x88$\xd6\x80\xb2\x89\xd4\x9e\xca\x93UV\nX;#\x82(\x9b\xff\xf7\xc3\xa0\x0cw\xafhl\xe7x\x83\x92\xec'\x8fJ\xa0\x12a \x95\x0em\x8a$'*\xf1\xc6\xb7co\x9e+\xcaV\xfa\xba\xe6>r\x8a\xa34\x83\x05\x0c\xa9\x0cI\xbc\xcaJ\xaeH\xf7\xe8\xa6\x89\xf1\xc2Z\xa8]/\x1f\x16b\x99\xe0\xb6\x94-\xd9J\xa2\xba\xcc\x9ct\x18\xe7-\xebF:\x89\x10\xa2-\x9b\x02\xa1d\xf9\\Y;\xdd'\xefwe\nY\xd2\x96\x91\xbbD\xe4Y\x94\x81u2e\xf2\x9e\x81\xdd\xe1+\x06&\x93\xdc\xdc\xe1e\xbe\x0f\xc0\x1fe5\x8d$\xe37\x94\x90\xb0\x0c\x0e\xb9\xbdn\x0c\x1eE\x8e\xcd\xfd\x8e\x96\x91\x9f\x04-a\x80\x96\xfd\xda\xb8\xb5\x94]\xb7\x0b\xbf\xcd\xfd\xad\x1a\xb8\x89,\x0cC\xc6\x98a\x96\xea4\xa5\x94Y\xebT\xbc\xf0\xa7]\xac4\xcb\x92%(z\xb4\xbb\xa1\x8a\xecMV\xa7\x0eP\xd9\xd5N\xd5\x8a\xe1%\xc9\x18\xa9\xc1\xc0}\xdbq1\xdf\xa8n\x15\x87\xe0\xc2\x05\xe79-\xb7\x17\xb5H\xadT\x13|\x92\xc8[\xe8)~\xed\":\xa9\xb1U\xa4\x9f\xd2\xcf\xce\xf1\x96j\x13~H05\xbd$\xea\xe3\xf0@\xb62\xbb\xc0\x86\xd5h\xfe\x16h\xc1\xec\x9e\xbc\xe3\xf3\x04\xe2\xa1l\xb84Ke\xcdvd4x7\xa9\x976\xecl\x1eM\xc0\xbd\x14hOB\x13\x8a\xe6\x0f\xd3\x10+p\xafc\xb7-e\x95\xb8!\xb7\t\x1a`\xa3\x04^t\xb0\x11\xca\xe8\xe3j#]$\xea\x99\xcb\xe8\xd3\x80\x8a\x9fd=\x9cK\xfe.\xe4\x8ap\xa1 s9\xb0R"

ref = b'BZh91AY&SY\xc1\xca\xefZ\x00\x00\x97\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xd0\x05^\xf2\xb7\xb7\xbbz7u=\xda\xee\xd8\xe9\xdd\xed\x94u4\x03\xd23I\xea\r=&4O$4\xda#F#\xcadd\x07\xa9\x93\xd2zF\x86hh\x9aha12\r4\xda\x8d6\xa6\x83G\x90\x8d\xa16\x8d\x11\xa6\x83\xd26\xa6OS@\xf5<Q\xe8\x8d\x92l(\xf5=\x08h6\xd4\x9e\xa6#1\x1aO\x13CQ\'\xa8\xd3&\x9bS\xd12z\x87\xa4\x1a\x01\xea\r\r\xa9\xea3SF\x13\x1aC\xd1\x94\xf51=\x13\xd4\x1a4\xc9\x91\xb12\x9a<\x90\xd9#jm&\xd4z\x83\xd4h\xf0\x89\xe53I\xa6\xc8\xc9=5<!\xa2c\n\x01\xe94\xda\x81\xb5\x0c\x9a\x18F\x80\xd0\xf5=A\xea$d\xf4\xd2hm&M4\xf52b=\r\t\x93\xd4\x18\x04\xc0&\'\xa2\x18\x98\r#F\xcad\xf4&&&d\x00\x00jbmF\xc8i\x18ji\x80\xd4\xf4&\x11\x82i\xa3j\x1e\x93\t\x84\xc3"1\xa0\x98\x03ML4\x0e\xa6\x8fP\xd0b4zM\x0fS\xd44f\xa6\x87\xa2by5<I\xe6\xa9\xe4\xd3L\x9a\x9e\xa6b\x9a{R~\xa16\x84\xf52=5\x0c\xd4\x1eSOP4\x0fS\xd3P\xc8=&jh4=F\x86\x8fQ\xa3Lj<\xa6\x9ad\xcdF\xd3P\xd3#\xd4\xf5\x06\x83CF\xd4d0\x9e\xa1U6\x9a\xa6j6\x93\xdaQ\xe4\xf5C\xd4\xd0\xda\x9e\x91\xea<\x89\xfaP\x1b\xd4\x83\xd4\xd0i\xea=F\xd2\x03M\x06\xf5F\x8d\x01\xa3j2d\xf5\x0c\x83M=M\x1azM4d\x1e\x91\xb5\x03\xd4\xd3\xf5OSM\x00\xda\x86\x9a4\x1azM4f\xa3\xf5COB\x1ay\x10hi\x90z\x8f\x14\x0c\x01\x84\x07c\xc5\x80X\xb2\xe7p\x05.\xea\xa9\xe4\x1f\xc3\xb4\x89m\xb9U\x7f\x16k3\xca\x92\xd1X\xc5lN\x832\x0f\xd9\x93T\x16\xa8N0\x1d\xa5\xc1\xce5\xa3\xc0\xe9\xcd$\x94\x08W\xd0E\xd8\n\xfa\x07\x0c\x8c4\xdd\x82\xa3Nd\xcb\xc5;\'z\x0f\xd6\x01\xd0E\x01\xdd\x89\x00\xd8\x8d2\xe6u\xc9L;0\xc0\xb0.c\x92V 51\xf8\xde\xc51\xea\x87A[\xeff\xba\x86=\xaap\xacz\x08\xaf{\xd2\xa7c\x8c\xc8U\n\xfb/\x94\x98\xc4\xa1q\xb4\xfc@FT\x9fK,\x9f\xfa\xfbA\x00\xf0\xda\xe3\x0f\xd0I\x81c\xc2\xbf\xe4\x87\x80Q\xe1\x92\x83\x96\x8e9"D\x87\x9d\xee\x7f\xe1.|U\xe8\x18\xaa\x90\x17\xd0\xa9K7\xf2t\xba\xdb\x99"n\xf0\xf1yw\xe8\x8a\xabK\xaa\xda]\xe0\x9e\x0c\xb0\xa0\xe7\x93\xc5\x83\x0eLY\x167\nh\xa7_%\xf9\xb1\x8f\x85#\xbf\xd4\xe1\xb2v$\xfc\xb9\x00_\x1a\xabY?\x14\xe1c\x07G\x15)\xd8$e\xd1\xca\xfd\xe0\x92\'\x0b\xa7\xbd\xc4+\x94\xddO xbS@\xce/\x84\x98i2\x99Hd,\x7fTI`\xce\xbd\xb5\x0c\xc9\xda\xdc-\xa0\xfb\x06\x1c\x06\xb1\t W\x8b\x0f^\x06,\x13\x03_\xa5=\xa1\x88z\xca\xdb\x06\xac\xa3N;-\xdfCR\xf9\xcc\xe2\xb36<\xc8\xce\x03\\q\xd3vvL\x90\xd9\xd8\x13i\xcf\xc7\xf9\xef\xf8\xaa\xdeW^W\x03l\xa2P\xba\xac*\x0bW%>r\xef\xe5\x96\xcdAr)9k`\xef\x1b\xcd-5e\xbb\xbaa\xd3z\xc0\x85\xfd<\xe5G\xc7\xc0mR\x1cBm\x91T\'\x0f\xcd]R+\x81\xc95\xf47\xa0\x89\xbdU\xd0I\xe9J(^\x1f\x1e\xc7s\xd4\x82lh\x0f-*9\x9eD{,WZ\x91=v\xce2\xc2\xfeL9%U\xac/\n\x87\x08\xc8\x0c\'o\xb4\xfd\xd9:\xc2Q\x08D\x07\xc9\xc1bff\x00A\x0c\x86g\xae\xe1\xe6n\xfa\x8dL\x89\x07\x87\xb76A\xcbU\xd3\xba\x85\t\x8c\xaf5\xb0j=2E\x89.\xe8\xa7\x802S\x85\x9d\x19BL\xae\x1e$+\xaf\x0f\x18 \xcd,\x9e\xfdv\x1b\x04\xdab\xf72e#\xda\xe4"\xc7\xec\x13\x0cG_\x12\xee\xa7\xefn>\xb1\xf8\xae\xed&\x82%x\xda\xb0\xc0\xbf\xaa\x9e\x8c\xb0k1\xcei]\xea\xda[\xa2u\x08\xf8\x8e\xd5\x95B#\x93\x0b\xe9\x1e\xac\x84J\xa3`I5\x01\x0e\xea\x9c\xa8"E\xabm_2G(6;\x8fvT\x0bS|\xaa\x90Y\xe8\xf7,\xf3W\xfe\x01\xd4p\xbe.\xe7o"\xae|\x16\xe6Q\xbe\x8c(\xc9\x8a\x03cc}v\x1a\xb9\x89\xeb1\xfc>p\xbc|(\xc4X5\r\xe5m\xa2\x84\x0c\x82\xc3\x1aX\xf4\xfft\x9b=n\x00\xdd\x9b\xf3\x81\xebU\xe6@\x9f/\xa1l\xc1d=Og^z*\xc5\xf6\xf5\x86\xddsR\xdfR;\xf4\xa2\xf8\xa3g\x19\xa7-\xb1\x1a\x02\xfc\'\xd1*\x82\xd9t\xe4f\xee\x9a\x92\xd7\'\xf8E\x9bo\x88\x1c\x04M\xf1\xd7\xc3\xea\xcd\x85\xc1\x958y\xa8\xa5t\xf0\xf0z\x04\x86\x19\x85\xcf\xa09\t?P\xff\xfcxZM\xd2&\x89\xfb]\x00~\xaaZq\xb2\xe9\x0cN\x94\x8bw\xd8\xc21n\x96\xbb\xc7\xc5\x8a\xb0Qu\xcbo\xbc\x03\x7f\xde+}\xf2P\xc2\xb7V\x8c\xc4\xd4\xc9\x89L\xf9\xa8\x93:\xc0\xdc\xc5\xb4\x170I\xca\x18-e}\xc9\x1c,]\x98+<\xd6L\x08H\xddB\xe41\x9c\x1a,\xec\x0b,e\x1f9\xa3\x84\xea\xc5\xa5\xf94%M|Ur\xe2|[b\xb7\x04\x1d\\\xda\xfa\x98(pI\x03p=\xc9|\xde\xd0`J\xe0\x03\x94;\xb7\xb6\x95C\x10\xe4\xd1\x94\xf0\xfa\xe2Z\xbc#\xdde\xaa_/Q\xd0\xa2\xb7\x89\x91\xbfo\xa5a6\x86\xf6\xa1\n\x992\xfe\x9dF\x08\xac\xd6W\x8aI\xaa\x1b:\xcaj\xde\x8a\xb3s\xa1\x1cd*\xdf\x9e\x16e!#\xa8\x0b\xe4\xe41\xa9\xefa\x17m\x93h\xd713J\xca*\xd1m_Y\xe9\xcb\xec&4\xa2F\xa3\xf2)\xd9\xd4\xc1p\x82\xf0\xe0\xd0\xeba\xd0\x88\xdc<\xf8VbQ\x96<g\x80$\xcc\x92\x0e\xa3\\\xc1#\xbc\x106<\x82\xdd\x0fU\xe8\x97\xb1\xa8\xd9Y\xf3[\xbd\x02\x1dh\xd1\x89\x92\xf7\x051\x87\xf6 <\xe8\xe7[\xdf<D3\t\xc4$\x00\xbc6"\xcd\xe9\x8e\x81i5X6p:\xb6w\xe7\xd6~\xce\xcbHK\xd2V]\x82bgK<\xdcV\x8d\xae\x9e\xc2"\x0fM\xccO\xc6\xd7\xc0\xeb\xd1\xbc\xed\xc4\xa5z\x9d\xb5\x10\xa7\x1aG\xe08V:\xcc4}o5\xcc\xdcf:\x1b"\x83\xb2b\x00\xed\xa3\xd5>\xa9g\x9f\xda\x89\x00\xc7\xef\xa8k%|\xc4[\xaa\xa1\x87\xacR(d\xc1\xa5\xa3>+B\x03\x8b\x02\x88K\xbf@\xc71y\x8c\xe1P\xb6\xca\xeb\xc3d\xeb0\xb8\xd7\xaaE_\x07\x8a\xa6\x9c0e^\x9e\x11\x89\x98\xadd\x16d\x17\xbe\xc8,:\xfbO\xc1n\xb3\x1f\xb8\xd7P\x14ws\xd2G\xd1\xd2\xa6\xc4\x0f.|i\xd0v\x85\xc9b\x831\x10\x0f\xe6\x12f\xe7\x8c\xe0J\x02\x08\x14lt;^\xe1\xe7\xc4m\xc6\xa9J^TP\x8bp\x83Q\xbe\xb1\xec\x03$DN\x89\x0e\x8ayc\x17\xfb[OSLrTt1(\xeb\xe5\n\x16\xbdDF<e\xdc\xed\xf5x\xe8/\x92\xd8A\xb4G\xea(\x0c\xe3#u\xb0\xb9_\xd5\xd4Vn\xf2\xd2U\xe0\x1c\xb5G\x90\xbf\xb9\xd1\x06k\xa4\x98Z\x10\xfa*\xb9\xb8mGHC\xcd@q"K\xa4r\x80=\x86\xd8Sa\xef\xcbO\x1d\x9b\x15\x1e\xad+gi>\xa4\xcdv_d\xb0K\xc3!\xa6&e,]\xd2\xe8\xc2\x0c\x83*\xbeB\xbd\x14\x15\xbax\xc9\xfe\xd1N\x1b\x8b\xc9\x8c[\x00\xd7\xfb\x12\x08\xc7\x99=\x06A\xa0\x1d\\\xcd\xc7<\xb9E`\xb1\xe9b\xcf\xa8\xa25\xee\xa9!\xd9\xe1ZY&<Y\xbb\xdf\xb7\xe0\xec\xbc\xd9\xdf\xa2\xa7\xaf4\xeb\x1ax\xd6\xa2(jo%Q3hS\xfc\x1a%\x01\xcd1E\xcb\xe9\x10\xe9O\x838xj\xafW\xe9\xdas\x1fO\xde\x9e\n\xc5l\xd3%\xa5\xf9\xe00\xe3/\xf6@\x84$K\x92\xd4\x0cM\x87\x8f\x96\xf9<\x9e9M=\n-\x08\xec<j\x9c\xf3\xdc5\x98\xe4)s\x0b\xeer\x81\xa4\xe9w\n\x03\xcb\xc4\'\x8e\x95\x02\xc8\xef{\x06F\xec\xdfv\xad\xfe\xdf\x00\x1c\xb5P-\xec\xb0\xb7\xd8\x7f"\xa2g\xa9\x01\xce\xc6\xd6\x94Y\x85fNG\xaa\xce\x80ui#C\x19\x83\xf2\xdd&\xb8g\x10\x15\xaep\xfc)\x11\x12\x95\x96b#K\x8f_\x0b\xfc;\x82[]\xc3\xcc\xbe\xac\xfe\xa8\x97\xd37hS\xd4\xc9`\xc4\xf1X\xa1\x18F\x0b\x81\xff\xa81\xe7i\x14j\xfe\xda\x04_x\xc4\xff;\xb6\x12\x99\xe3}\x10\xd7\x1b\xf3L\xcc\xce/\x077rs\x1d\xc7\xf9mI\x19Hb\x9a~e\xe0\x95\xa9_\x03\xbbX\x8f\x8c\x84%\xc0\xb2;\x91}\x08U\xd8r2\xc6\x99$H\x15c!\x8a+D\xb1-\x8e\x83\x84V6>\x1a&\xb2BXJ\xe0\x9b5\x9a\xe6*\xa61y\xb6\x80{\xfe\xb2h\xbf \x85E\xca\xd4\xa2\xc9i\xf6.\xd8\xc0\x93\xb0\x8f4(r\xa09+$":/\xda\xf6\xaa\xc4D\xa2\xde/e\xccr\xf3\x1dE\x82\x0f\xda\\\x86\x03+sg*\x9a\xe6O\xd2=\x98x\tqB\xcf\x97\xa8\xe8\xf8}*Wm\xe9L\x03!U\xf9\x029>\xa4Lm\xbb\xda\xdd=j\x97\xd5]\xd5|4q}\xe6\xef\x169Y\xad4\xa0\x82\xc1\xeb\x80&\x16\xf4\xd3\x9f\xd1]\x84\x99\x83}\xb5\x18\x13\x98\x18\x9da(\x1c\x9b!\x94U\x11\xa0\xb7d\xff\xe6&B\x0e\xc6=\xfd\xba\xda\x195\x18&\xda\x97a\xed\xfa.]/\xa3\xe0\xb5\xc78\xca,\xcd3\xaa\xaf\xe8Z\x04\xd7\x95\x92\x98\xb8\xc2O\xf6\x9d\x8bP6\xeb\xaaC\xda\xdd\x8d\xcb\x1c\xd5pw}\x1cL\x003\x8e\x95yXJ\x03"fj\xe1\xe4\xee7\xf4\xa8\x8b$+%N\xa3\x9a\n\x02h"\xc6\xa7\x94\t\x14p\x17N\xe3\xbb\xbc\xb8\x82u\xf7\t\x92\xc0 \xf9\xff\x17rE8P\x90\xc1\xca\xefZ'

all_cram = b"BZh91AY&SY\x8d\xea\x1dy\x00\x06\x84\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xe0\x0e\xde\xf7\xc3C\xa7\xd0\xfb\xdfK\xbb\xdc\xf6_=\xde\xdd\xb1\x9f\x1d\xeb\xdc\x1f'\xde\xb7m{\xd7\xce\xbe\xf6\x1d\xad\xed\xad[$\xb6\xbe,<\xe1T\x13L\xa7\x93Q\x99\x0c\xa9\xe6\xa2y5==S\xcc\xa9\xb2i\xa6)\xe0\xa3\xf2!\xa6&G\xaaz42f\x93\x13\nzcR{$\xc2\x9e)\xfa\xa3d\xda\x8d3S)\xed&\x862\x9a4\xf512z\x8c\x8d\xa2\x13\xd3\xd0\x13\xca\x9e\xa7\x88f\x83M#L5\x0cI\x93\xcc\xa6S\xca\nH4\x14\xf3U6hF4h)\xa6\xd1\xa9\x9a\x9eSjz\x1e\x8d)\xe6S\xcajy\xa9\xea\x87\xa9\xe5\x0f\xc9\xa5\x1ehD\xd3z\xa3\xf4\xd4L\xcd\r&\xd2Oh\xa1\xec\xa7\xa2\x9f\xa3L\xa9\xfa&\x98\xa3\xda\x14\xfdC\xd5=\xea(\xf5=G\x854\xf5<\xa7\xa8\xdaOA\xeaG\xa2f\xa7\xa8\xc2{T\xf4M\x07\x91\x92hT\x99&\xd4\x99\x95?S\x1a'\xa1\xa0M0\x8fA6\x99MO\x1aQ\xeal\xd0\x9b@\n?\x15\x1e\xd2f\x83&\xa6\xa6\xc9\x80\xd3H\xcaz5<\x93zS\xd4\xd3\xf5M\xea\x9e\xd56S\xc1OS2<\xa17\xa6\xa9\xe5\x0fI\xfaSj<#CMM\xa4\xfdL\x9a\x99\xe9\xaa~\xa8\xd0\xf4\x8fS\xd1\x94~\x94\xf2\x82$\xa6\x12oT\xd8CL\x89\xeam\xaazmD\xf5=\x01\x91\x92z\x98\xc9=\x14\xf2\x8d\xe5O'\x91\x1a\x8c\x0c\x93'\xa6\xa3\xf5\x03I\xe9\x18\xd2m4\xd0=LCOSS\xd1\x8d!\xa2mOM=)\xea\x1b'\xa5?H#\xd3!\xa3D\xcd=Bdzjm5=\x02i\xfa\x98\x83M\x01\x08\x82\x14\xf6\xa1\x89<\x93\x135\x0cC\x13&i\x1af\xa7\xa814d\xf5\x00\xd0\x1a14mA\xa3 \xd1\x88\xd0h\xd3i\x00\x1a\x1a\x1e\xa6 \r\x1a\x00h\x06#\xd4\x1a4\xd3@\x00h\xd3 h\x1a\x03@\xd3@bD\x80\x8c\x93F\x9a2h\x98\x9eI\x9ah\x08\xd0\xc9\xb4&\x13C\x08\xc4x\x89\xa7\xa4b\x1e\x88\x1a\x18\x9a`\x9aa0\xd1=\x13#\t\x84bb\x1a\x1ad4\xd1\x90b\x1a\r\x18\x9a\x03CA\xa3 \xd3OBzG\xa8\x17r$\x11\x10\x00>9=|\x89\xe0\x00 \r\x08\x15\x02\x03[\xd5\x89\xadz\x85\x07$\x128\xe3\x92C*Jd\xa3\x04\x89Z^\x9f\x89\xc8-\xa78\xedb|\xe3\xf2\x02\x81\x1e\xfb\x91\xf2\xbd\xc2cG\x19\x15\x94\x8a\x1c\xb6<\xd4\xbed\x88\\\xa1\xef\x89Q\xb87\xaf\xad\xfc\xac\xd2\"\x80\x04\x0b2 \x11-Y\x91kx\x15\x95`DY%C\x07\"\x12!\xe8\x80\xc8\x0e\xd1\x16D\x08\x91XP!\xda\xa0\x08\xaf(J\xcf\x11\xc3Uv\xda:\x10\x00\x80\x81+/#kX\x96\x80\\\xfe}\xf4P\xf1.\x02!\x15\xff\x06(\x06QP\xdc8T\xe4\x1a\xcd(\xca+s\xf4:\xb7\x134\x8b\xa0U\x86\xef\xab\xf0Y+\xac\"%\xfb7\xdb\x8c\x02\xde\x80uB\xcc\xa7\x87\xa5(h\xd0p\xae\xb1\x10w\xb3\x0e\x06\x1c\"\x1e`\x1e\xbc\x8e\xb3G\x18sTSl\xaeP\xa7\x84Q\x04\x1eRI\xbfA,X\x89\xed\xaa@\x86\x10U\x16u)p\xb2E\x16F\xb0\x84\x99\x1b\x82\xb6b\x98\xccLp\x9c&\xd0\x1e\x9d\xf1\xb1\xca\xf4\xb4\xd4\x9b9/\x8a\xac*\xe0\xfe\x9d\xd2-0;\x10\x9f4\xc3\x0f\x1b\xd8q\xd5\xef6\x0e\x9dC\x12\xa2\xa0\x04\x98X\xb7\x0f\x9eL\xda\xf6\xb1rk\x02\x1crDr\xc6\x9e\xd1>x\xdc\x08P\xa1\x9aB\xe7;\x1c\xac\xba<\xd4\x07&\xa0$\xeb$\x92p:+Q\xb4\x83m\t\x99\x8d)\x9a9R\xf3\x01yu\x11\xae\x9d\xdcr$\x12\x00\x8a\xba\x94\x01\x03\x92<\x1c\xbf\xdd\x9b\xa4\xf2~\t\x99\xfc~\xd3]\xad\xec\xafs\x1c\x19=v\x86R/\xd8\xa1\xeb\x9d|\xf4\xf1sG&\xba$\xa8!\t\x8a\x8fk\tb\xf6\xb8\x9du\xf0\xd6\xc1A\x17v\xc3\x03KLZ>\xa4\xc2\xf8H\xfcv\x19\xff\xdd%\xb7\x1a!\xc0\x9a\x0er\xaf\x14tvs\xb7@\xf4(~\x8a2\x00\x8d\"\x04y\x17Su<\x8f\xab3\xa70S+\xd8\xed\xc7\x85j\xcf\xb8b*\xe2p2\xd8;\xd5^\xa3 t\x1em\rnA\xb1\xa0k\xc5\xa3\xcd\xdc\xf2\x95>\x83\xc8\xb1M\xb1 \x0f\x85\xaa.\xee\x00\xc4\xf4\x1a\x9f\x80v9M\x89\x00MC\xc6\x94\xfe\\\xa2E\x98\xc3;Ee\xbcX\xe5\x83\x9dC>L \x80\x1a0\x08 B\x07\xc7\t\x18\x04\x16\xe1\x8c\x8de^\x81\xf5\x99\xc7\xc9\xee\x1bX\x9e\x86\xe12\xf5\xec\x9bg\xea`5\xc7\x06@\x94\x80A\x04\x04\xad\xb1\xc1\xc9\r\xaae\x11\x06)m\xe6\xdd\xdf\xa2=\\\xc7w\xcc\xc4\xd4\x8d%.\xc9\x9b\\\xda\x02A\x96\xb9:\x80\x02\x06\xa3-\xaf\xc1\x17\x94g\x89\x13\xc6\xb8^\x00\x03\xce\xb9\xa62\x93\xfc7q\xf4z\xad\xee\xc8w;\xca,\xccx\xc0\x80k\xee\xb4H\x17\xbeO\x17\xfd\xb8EY\x1eh\x0c03g^A\xa2 \xe5u\xba|5\x81\x1f~\x96vou#\xd1\x8d\x9e\x9c\xb7\x80\x14+fj\xefg\x8e\xc4\x18#\xb2\r%B;\xf2c\x16\xf2y\x144\xb5\x90zp\x1a\x15H_\x1fy\xbf/\x8b\xbc\xec\xe7m\n\xebl\xfa\xd1\xd4\xcb\xfe&dx\xaab\x85\x14l_n\x0b\xa9}\x0e\xe8\x88\x80\x06\xeeng>\xfa5\xab\x97\x8c\xb7~\xe8\x0f\x9fJ\x0b\x030Tg\xc0\x87\xbaz\x12\xc0\rV\xd7\xa9\xdf\x07\x19\xb4\x86\xe5\xc29P\x00\x94\x18@\xdb\xf7H\xbbXc\x81f\xffSH\xadn\xa9\x1b\xfa\xc7\xc4\x9f\x8d7\xb7\x1d\xe2\xe1\x10\xb1g\xf0\xc6\x1e\xce\x9a\t\xb9Y\xcfCQ\xa7\xa14\xe54\xb5\xf0\x8c3\xb9gu2\x1e\x94\x05\xfc\x08\x02$\x986r\xa5\xd0\xcc\t\r\x08\x1d\x18\x98\x7f\x9e&2a\xcc\xf2\xb4\xb8I\x0cc\x104\xd1!\xe2^y\x1d\xde\xf5o\xd4\xfan\xac\x1cEZ\x1d\x16M\x03\xd7#\xc4K\xe6\xb0)Q}\xd2KU\x8e\x03\xb4\xc3\x0f\xbc\xd0\xae\x10\xa1\x93\xf8\xfbi\xe3 \t\x03B\xcda9\x951\xd7u\xd9\xe6b\x9b\xf0\xe1\xcc)\xc1\xf8J\x0bw\xe6\xae\xa9\xdb\x8a`\xde9D\xe8\xe5K\xb8U\x13\xf0\x85=\xcbv\xe1x\xf0\xa7\xb3X\r,\xc8\xdc\xac\xcc\xb1\xac\x1a\x8cm\xff]A\x06\x08C\xbc\xaf\x97%\xb1\xb7C\x14\x1dX\xd6\x0eA\r\xcem\xc0[\x96R\x12Kr\xb2?\x94*!\xd87\x97\xd5r\xd46\xdc\x0b\x93.f\xa7|\x0ejj\x1c\x81\xa1\xbf\x83\x03?\x93{\xb7\x89:0\xab\x85\xc8\xf0A\xd5\xa8%(\x18n\xe9\xd3\xd9\x16\xdb\xe0\x0c\x02\xa6E\x84\x03:\x86\xe1+\x94\xecr\x8a\x16&\xe4t\xa2\xa3\x82\x9d\xcb/S\xd5Y\xd0G\t\xe7\xa5@\xca\x9e\x86I+\xf1/\x9b\x00\xc6R\x9f\xac\xcf\xb1x\xf1\xfd\x00\xcbtUq\x829'a/\x00\xb5\xc8\xc7K2xy\x126\xebk\xf7\xdai\x9f\xda\x80\x8a\xf43{\x88\x82x\xe7\x04t4\xc6\xae\xca=\x8b@\xfa\x82Ey]z\x8bJ\xd30\x0e\xf8\xf1&t\x0e\xd1\xb3b\\r\xe6@\x01\x1c\x19E!\x06DJ\xde=\xa2\xa1\x9e@Z\x98\x84\xa3\xaa\xa8h_\xd4\xa4\xd5es}\x0b\xee\x1d\x9d\x91J\xba\xca\x93\r\x16El\x17Wj\xbfJ\xddB\x1e\x17\x88\x9e\x89\xf8\xc6\xebHo\xcc\xcf\x8e\xd3tJ\x17%d`\xdf?c\x84\xe7<\xd8x\xd4\x84V\xd6\xf5\xad\x86\xb9\xd8\xf8\xf9\x9a\x12\xf1\x0f\xea\xc6\x00\x90O\xb0\x01vf\xdc\xe4\xbd\x8b\xe9\xe9\xdd\xe7\"`\t\xbc\x96O\x83h\x16\x06\xa9l\xb0\x8c\xa4\x8e7\xf0\xee\x95\xfd\xd3\xe3\xe5\xd4vL|n\xa4<\x8d\x91\xf0\xc7d\x9d\x1a\xce\x80[4\xc6\xa1\xe1\x14\x11\x0eD\x05\x9f/@\xa8\x8c|\xf0\x82{\x10&E.\x96\x80\xcf\xaa\xc9Y\xd9]V\xf1q\xae\xe4i\x05\x17\x081 \x9d\xd7\xfb\xcf\xd6\x04H#\xa4l\xfalw\xe3mHJ\xfeK\xd8\xa5\xbe\xd15\x15\xed\x07\xce\xe1\x93a\xbb\xc7\xe6,\x8cy\xb9\xc8\xc8\xe8\xc8\xd2\x1b\xc9N<_\xc7\x83\xcb\xd0okf\xc2\x0e\xe7\x8bY\x06\xec\x8e]\xe64\xa0<!\xd5Oc<=\xaa\x8b=\xcb\xb5M\x0c( ^{R\x96:\xa0lj\xf3)\xe9Y\xa6\xb6\xed:\xee\x9e~BK\x02Ed\xde\xf6\x8fc] F\xf3\x98\x83\xf7\xa7\xf6\xf5\xdc\xa4\x87f\x1a\xfe\xca\xf1\xf4\x0b\xe2\xab\x85\xb3}\xa7Q\xf9\xf1'\x90e\xe73\xb9^|\x8c\xf5f\xba\x05\x18\xa9\x91%\x9e\xe3\xd2\xc1\x16\x04pQ\x08w\xd8\xec\xf9\xe2&d\xfd\xf5d\xb4aV\xa2[\xcf\xcd7aTX@\x18W\xa4\xad2\xd68\xa4S\x8e\xe1\xa6M\xf2Q8\xd0\xe2\xb3\x9aX\x84'\xd9\xa2?u'!7ug\xf2\xac.[\x8e2\xe8\xdc2\x1a\xdf6\xda\xc3M\xe9`\xb8>\x14{r|\x04\x93\x836\x8a\x11\xbf\xd0\xb4\\W:\xb0\x90\xc9?\x1d\xa3}\"X\x0b\xcb\x06\xa8A\x04\xa0E\x02\xdfaw\x8d\xc8\xfapK\xd8\x943\xb3\xa7^\x8fY\xde\xa0Al\xd81i1\xc4\xdc\xa2\x05\x18\xc9\xa8+\x1e\xc1\xa4\r\x8cu\xa5\x13\x81'\x88:1Xj\x1b\xb8\xfb\xb8\x9e\xa0\x90_0[{*$\x8dQC\x92\x0c)>\xdd\xcdt\xe1|S\xb6\x8a\x9e\xdd\xc3/\x88e\xc7K\xaa\xd1\x94*=\xfe\xa1\xe4\xa8D_\xc8{=K\x94h\"\xc2\xd8\x865\xd8\x9bIN,\td\x84\x0b\xe1\xcb\xf9\xb7:Z\xdaR\xbba\xcb\x900%\x1d\xfe\xa2\x1c\xfb\xc3E,\x1b\xf2F\xa2\xb2\x91=\xebu\xa4\x1d\xdc\xa8\xb4\xba7\x9e\xbbU\x15+\x1d=\xd7\xa7\xe8\xb5\xf3#\x07\xb2\x11\xfeK\xa3Y\xf0h\xe3\xdaHv\x1f\xa2\xdf\xb3\xb9\xe8\xb4\xb7c\x0bu\xcf\x01\x08\xfa\xb00\x9d\x83\\\xf9\x89\xe9\xd0\xb9\x8d\x82\x1b\xc3\x0b\xab\xac\x8bH;N'K?\xef\x9bZ\x9a\xf8~\x82\x1b;\xc59\xcb\x0e1\x8b0\xd5\x1b)\x8b%'6\xfb5\x0fPV\xb6\xdft\xfeR\xcb\x8eo\x89\xf4\xb2G\"\xd4\x91\xeb\xc1\"\x8d\xbfA\x99\x1d\x8bs\xfa\x89\xdc\xefm\xad-\xb3\xaeO\x96\xc3\xafq\xc1\xb3ufX\xc5M+U\x93\x87\xd5giu\xc1\x1f\xaa\xe9\xadC\x81\xe4_+\x83[\x1c\x1b\x91\xf7\xfa\xb3XNh\x05\x186w\xb6\x17oB\xc9u\x0eH\x99{\xf0\xeb\x1ebc\xd6\xd0\xc5Sl\x98\xd1=\"\x00P\xa9\xc1\xa9\x80a\xd4\xe2!$\x0b\xc1\xb0\x80\x0e\xfe.\xa3\x03\xe2\xdfA\x14%\x1c\xdc&\xe6\x8fLX]\x0f\xe5\xdbu\x19:\xda\x84\r'\x89\xc04\xa6\xb4p\x06\xe6\x18\x14\x19\x92Hm\xd5\xb2\x99\xb2\xe6\x02\x0c\"\xe2Z\x06\x00\xed\x11 ]\xe0$\x96\xfb<\n\xae\xb3s\x00\x8c\x82!\x08i\"\x82\x19\xbb~\x8e\x1c\x9d\xbf\xf3\x8c\x10\x82\x8eyQ`'?M%E'\x8a\xe6\xde*\xb2Q=\xd9\xc7\x8a\xc9\x9d\x16\xa7\xf5\xe4 g8\xe0\xe6V\x9a\x05\xcaUB\x91\x119Q\x0c\x0fLy\xbc\x1f\x94\x85\xf2\xb0\x87\xa5\xea\x04\xa4+\xe3ZS\x88a$B\x9b\x91\xa0\xe5\xf6\xfc(\xe9\xf3\xdc\xd9\xdegh\x12\x89\xd4\x1f\xb8\xf5\x94\x9a\xb5\nYo\x87i\x9c\x08\xd4\x88\xf5t9\x01ig\xfe\xa3\n\x92\x8a-BMi\x9b\xc2n\xe3R\x15\x8dg\xb3\x9b\xde@\xc7\xf7,\xbc\n+j\xc2\xff\x99\xba\xab'h\xde\xddF\x9e\x12\xaa\xfd\xd5v\xfbj)\x16\xbb\x0c\xb6\x99Q\x81\x183\xf5Z\x7f\xcb\xae\xdd\xde\x83T\x1e<\x13\x8f\xc6\xd6r'\xcb6\x13\xf3\xa1\x0f\xb2\x10\x17\rj\x0c0\xe4\x14C\xc0@1@8\xc5\xbfL\x17\x83\x040\x82I \x83\x03\x98\x9f\x16-\x84\x7fw$\xb9\xecd\x9e\xc4\xc2\xb1\x83\x969-W\xa5\xeb\xc4\x14R\xe06A0\x13,\x07\xf6\x87\xb7\xbd\xe4}\x1a\xdd\x86\xf1cdh5\x12,\xb2\xd1)*/Z\xcc\x15\xa0\xceW\xcb\x06\x0b\x8d\x84\xdc\xdbg8z\xff~\x97\xc7\xceZ\xf9\xd3\x9e\xaf\xd3\xf4\x0f$\x0bZ\x07(yTUt|J\x05M\xa9\xbc\xdf\xba\xbd\x8df \xd8\xdfB<SF^t=\xed\xb6\xc9\xc0\xdcSo\xbd\xbc\xd7\xd3`\xf1\x8e\x0e\xdcc\x96\x1f\xa8\x185\x7f\xb9K\xd0\x91]\x00\xa2m\xc7\x96fJ,Q \"@{\x02\xc0\xf1\x1eb\xab\x04\xb1\x82\x03\x16\x82\xd8S@`\xee4y\xb5u\xf2\x83x\xa4\x19\x94\xe1\xd9\x9e\x9f/\x07Q\xc5\x91G\xe3mX\xed\xb1\xe4\x87\x92\x88\xf8\xe6m\xcd\x8c)\x15\x1aI\x8a\xc3\x8c}\x7f\x80\x97\x92|\xad\xe5\xb4i(\x1c<T\xed9\x0b\xe7\xf7\xca\xb4\x05h\x88\x01\x00m\xf1\x13;\xdaS\xe4\xd8mi\xb3\x9b\x84\xdf\x1f+\xef\xc9Y\xc02\x87\x0c/&\xe2\x1eo\x8fps@\xcb\xfd\x84\x1e\x1c\n`n\x90i\x9e\xd0@\x044\x11\x0c\tq<W\xd55U\x86\xfd\xd7\xea\xea\xa2\x08\xa2Q%\x1cJ$\x983\xd2\x1f\x1f[W\x0e\x1c\x17\x1d,\xd9\xb2\xe8\x08(n\xdc\xf6\xe1\x1a\xf1\x0e]\xb8\xa6\x80a\xc6`\xe8\xf6\xb4[z\x12\xf8\xd5\xa7\xbf\xd4\xb7\x84\xb2\xb2\xad\xdc*\xf8~;F\xc7\xcdU~Ff/\xef\xbd\x9ea\x94u\xc1\xdeG\x8d\x17^x\xda\xdb\xd5r\xe9\xbbC\x8e;\x8f\x9e\x88\n%,!\x86\x18\x0c^\x8d\xf9\xf554\x0f\xd6$\xa2\x1cSJ\xcd\xaf\x9f=\x10\x9e/\x04<#a\x06(X\xcf\xba\xc6\xcfO:\x82\xb9\x8a+\xe6or\xd6w\xe3\xfb\x1a\x18\xd0\x02\xa0\n'\xddM\xfcT\xe3\x82 \x19\xce\xe1\x85\x8534\xc0\x81\xf2\xa2\x0b\x88\x1a$\x04\xcfn\xd53\xcdg\xf4g\xbf\xebL\x8d\x9evvi\xf0\xe6\x9d\xc7\xb8\xce\xa7J\x98\x82i\x1a\xe0\x05\"\xa0\x08P\x07\x8f[j^%4z\x87`\x9c\x9b\xa3\xb4\x9d<\x9c\xb62\x01J\xcb\x1b5\xa9]\xa2\xf6\xb7\xfa8\x16\xd0lV\x18p\xb9o\x88r\x1e\x80\x10 !\xc8\x00[*p\xc1\x88\xa7\x97\xd9Lj\x18<\xe2\x90\xd7\x00\x07\xef}\xe5\xdd\xef\xf0\xbaG\xe1\x89\x8f\xfd*\xb8J]\xde8r\xae\xddn&\xbaN\xc7\x7f\xd7\xde\x08\xbc \x12\x94\x05\x89@\xc6\xe9s\x98O\xb6\xcb^\xbbC\xf6\xbd\x9d\xed\x0e&\x9a\xd3;\xf8\xc5\xc4o6\nf\xd4\xe5L\xddM8z\xa8g\x1d|mR\x05\xf6\xbf\x9c=\x95R\x97d\xb6\x1a\x9a\xbc\xa8\xe6\xbb\xac\xf3B\x00\x06\xf0\xb0\x1ba\xff5I\xeaX\x0e\x13]\xc1\xdf\xd5\xa5K\xeex+n\xab\t\x95\xcc\x93-\xef\x92\x1d\x86+\xfc\xa0\x03\xb3\x08\x02 \x034flMa\xe2%}\xaa^4\xac\xb5m\xe9\xba\x02\xec\xa8\xa5\xc8\x88WQD\xb2\xae\xd2K\xca\x93=\x89o\xcd\x86\xd22\xac\xd15\x98\x06\xdd\t\xebW\xa9\x844\x0c0\xe8\xef\xa5\x92m\xaf\x0b\xc5P\xaaR!)\x10\x83\r4vKm\x0e\x8fQ\x87\xc2\x10\x00$l\x19\xc6\xe3\x1a\x88\xb2\xc0vh\x10l\xe3r\rK\x00!\x17\x93\xb9z\xb5\xba\xf9\x90\x8f\x85\xbavH\xbd\x7f\xf4P\xb6\xee\x0e\xbd\xbd=%\xbc\xe3\xf8\xc8j\x08&\xb0v\xd1\xd7\x19\xe3}6\x19\x08\xa6\x13\xb1bA\xd5\xc4\xfa\x8a\xc9<\x8a\xbe\x18\x9c\x16\xaf-)!!Ci\xbd+\xfcQ\xa1\x8fq6\xa4zG\x15<\x01\xeb\xc8\xa8\xb8\xb3x_U\xddv\xc2\xd4l\xda0\x1f\xd9\xac~O\xe448\xa8\xce\xb1\xe8s\xb9{3\x8bc)V\x17)N\xed\xbf\xb2\xf2\x99\xf7\xcbt\xde\xfb\x99\xba\xe5\x0f\xf9\xa76\x87U\xbf-@\"\x8f\xa7\x00A\x05\xe5\xbf\xa6a\x16~ST\xe3\xcf\x13\xd4\t^\x05\xa4\x94\xacT\xf10-\xcf\r\x18\x154f\x8c\xa4\x1d\x8e%\x16\xfc\xf6\xda\x80\x80p4\xda^\x9d\xf6\x04\x08:\x82\xae|t\xb3m\xa7\xc6;LR\xf6\xf1\xd3|\x15G6f\x989\x1an\"u5t\xcb\x8a'$k\xe1\xd1\xf32\xb7=\xa5\xe6\xf4\xdc\x828\xd2\x84\xb9\x1eX\xac\x935\xce\xd0\xcd[\xd0!\x14\x8a `s9\x05K\xb5\x83\"\xf7\x85k\xd9^\xeay2h\x1dm9\xd8|\x9b\x89\xfe\xb4\xef{\xfci\x81c\x85\xb6\xe2\"\x967I+l\xc3z\x8a/\x94bu\xf2\xe7\xcd\x07\xa0\xee\x8a\x92\n\xe6t\xb8\xbd\xf3\x15\x00~i\x16D\xad\x8e\xbeG\x83\xf6|\x9b\xaeL\xda\xd4N\xe3\xcd\tR\xda\x0e\xd1\xfe\\]<\x91H\xdc\x1c)\x1d^\xf5N48A\x04\x00\x03\xaez2\xf9p\xbd0%h\x1bL\x1d\x12\xe0\xe5\xce`f<\xd2\xdd\x83Yy\x1e\x89\x87I\x85(p \xf4_\xeb\"\xdd\xa5@\xd3.\x9b\xe1~\xef\xbah\x97\xe2{\xda\xd4@\x1a\x96\xaeLa\x0f\xd7\x90\"\x01A\xc9*\xb1\xbb\xd7\xff\xaa\xb0 \x05\xee+M\xd2\xc1%\r\xc47\x8b\x90Y\xbe&O\xa6#\xde\x97\x1a\xed\"\xa3\xfd.\x03\xaat\xc3\xf4+t\xf1aa|\xae\x89v\xb4\xefOCA\xf7E\xcb\xb6\x0f\x9d\xb13C\x03KZ9\xb5<\xd5\xc7\x96 V\xf9\n\xed,\xbd\xf9\xc7.\xc5\xf3\r\x85O]2}\xb7\xd7\x12BK\xbfE\x98\xb8\xe4\x19y\xbaL\xdb\xf7\xb9_\x16\x157\x81\xbby\xe8\xa9\xe0\t\xe40\xaf\x081\xb9\xf9\xd4\x0c\xec-|\xcf\xb7\xffJ|\x08\xc9NS\xb8X\xe1t\xf3\x10q\x04->\xd3\x9f2fg\x8c\x80\xa7u\xd1\xd2Z&U\xf3UQ\xbf\xbf\xe1\xec\x01\xfc\x86\xf1\x87\xaf\x1d\xad\xe8\xf9:/\x919P\xde\xdb\xf0\x18\xdb\xd9\x85\x8c\x19\xcc\xe8\xe1\x98\xb5\xc5#\xe0\x0c\xddq\x1b\xe2\x0f)_)\xc7!\x11\xde\xd6\xe3{\xa4\xb7\x10\x0f\x079\xc5<i:\x06\xde\xba\xed\xf0\x85\x94fZ\x91\xa8\x97;\xf9\x01\xdej\x92\xb88\xb8\xdd.\xf2\xf6\xc3\xd7\xef\x8f\xdb\x90[bh5\x88\x80\x1a\xe7\xbex\x1eeM\xdf\x05NC#\x8fL\x1a\xc2\xa6\xe8e\x98z\xdd\xc4'\xe5\xe4\xb4\xd2d\xd1\x8c\x99\xach\xf5T\xdc8Ye\xf5\xca\xe3?T\xa0\xfb]\x1bU\xdf~1\x19E\r4\xff\x8cGw:,\x96{8\x99\xd6\xd5c\xbb\x92\xd4<]QJ\x04\x8eC\xeeL{\x08}|\\\xd1A\x00s\xae\x93]\xed\x16gE.U\xb3[\x92\xf6\x0ez\xbb\xa3X}\x02>\xfd$fG\xb9\xed\xa7\xde\xda\xc1\x1a\x99@\xf7\x08p\xcd\xcd+\x17\xd0\xe3\xbet\x83c+\x0b\xc8\xd8&\x95\x94w\xc5\x8e\x91_\xc6e\xb0.||!c\x83\x88\xcf\xaaR\xa2\x00\x12\xbd\xd0,\x94\x85\x87\xdd\x00\xa8\x82[<$4\x1f\x8e\xc5T\xf5\xc4\xb8v/H\xc1\xb8\n\xc7 \x9f\xab\x8c}Dg\xd8u\xbfc\x9b\xfe\xe0\xcd\xdf8\xfa\x8f2\x10\x10\xf9\xa5\xfau\xec\x18=\x95\xaa\xee\x8bR\x81\xe2\x9d`d\xe2|\x84\x95\x111\t\xf2\xb0\x151\x8df|\x0b\xad\xb5b\xf26\xf3\x8fEfs^\xe0f\x9b\xf8\xc2\x00\x8bz\x03J\xb6\xf9L\x9aVM\n\xd7\xb4\xe6@\xfe\xe0\xb8\xfa\xe0\xf6\xf3\xe7\xf1\x8c0W\xb8\xa3\xe2\x92\x84\xc5\x1edZn>\x10\x8d\xb8\xf3Xr\x06@\xac\xbb$X+\xe6S\xdc,\\\xc3\x19\xbeg9\x1c\xeayVo\xbf\xf6\x19I{\x98\xa9&\x0f`\x1d\xe2e\xd1\xeb\x10\x03d\xec\xe7\xbam[6\xa9#+\xa8\xec\xed\x9f`h\xd8\xb2\xf9\x91\xa2=\xb6\x99&\xaa\xfc\xea\xd7\xddG\x1f=\xec\xb6\xef\x8a\xf7o\xeb\n\xd2\xd8\xe5\xcaxH]\x96\xc6\xed\xd8\x9b\xa6\xa7\x11\xbf\xcd\xa0\xb2\x9eth\t\xe8\xf3R\xc5\xce\xea\xea&\xc0\xcf\x0e\xd8\x1aZ.Q\x07K\xb5\xf5\x1du\xc7\xeb\xa8#\x02\x85\x8e\x91\xe2\xa6\xc6\np\xd8\xcaSW\xa1\xb7F\xb9\x0f\xb5\xa6\xd2\xcc\xfc\x85\x10\x12\xe3\xf0\xd4\x10#%\x9e{7Q\xaa\xb3\xcf\xca\xd8\"\xb6\xe1\xe2#I!P\xe0\xc4\xdc\xf8\xfe\xe3\x04\x11\xf4\xde\x04/*\xd9\xa7\x1c\xd6\x0b\xaa\xca1\x8c,V\xf1@(\xdd^\xaa\xb7\x9a\xf9H\x8bB\xf9[\x19\xef\x0fP\xf7\x8b\xb7+\xb3f?\xd3\x84o\xcd\xb9E\x0c(\xf6\xdcpC43\xa9\xb9\xb90\x07\xd0\xa4\x15\xd9\x82\xad6\\\x99\xd8 8\x01\xdc\xf6S\xc16\x8cx\x8e\xd1\xbb\x84\x81\xa5:\xd2\xa5*-\xcc\xf1\xf9\x86yl\xb8>\xcc%X\xc5k_zM\x1a'\x11\xe8k\x87{\x1f\xb7\x94(.\xe6\xf2\xce\xb0\x14\x95\xa7\xd2\xc9\xbbC\x00\x95\xd0\xa5)\x01\x81WF\x8f/\xbe\xd5L\xef;\x99\xd6\xde\x0fF\xd3\xaaT=\xa8\x13\x0e\xd6P)\xaaj\xe5\xe9!<\xd2\xc00{q\xc6,\xe6\xcf\xd8\xce\xb4\x83>\x83\x95w\x0bu\x01\xcc\x00\xdb\xbd\xff\xde\xeb\xd8\xfc0\xedl*\xe6\x96%q\xe6fX\xe7\xe2m\xd0\xc8z\xb8\xc2\xa2\xa0\x7f++\xb3\xec\x8d3\xb3\x1c_$~\x08\x13\x9eq\x9b\xaf\x1c\xb5\xea\xbe\xfb\xc2 >\xb4\x83\x1d\xea\xeb\x1f@\xfd{#\xf1\xc5\x05]\xd5yb\xd2\x92\x82?k|\xe71\xdb\xe2\xce\x97\xa2\xf0\xedE\xfe<\xca\xdc\t\xf2]\x9e\xae\x0b\x8cL\x17\x06\xfa\xabq\xa7\xda\xda\x15\x02r|\xff9\xc3\x1c\x9e\x81\xdb\x86\x85\xa3Q\x9aB\xff\xaf\x12/C\x7f\x89wxv\x90\xdc\x85\xc2g\x14V\xdb'd'3\xd5,r\x07H;G7\x90\xff\x1c\x01!\x98\xde\xb0\xd2\xdcLk\xad\xc8\x984\xda\xa0yO\x12\x1b\xca=\xf8\x7fz8\x9a\x8f\x8a\xea\xb8\x9c\x1cY\xe6q\xa4\n}he\x85?\x00\t$\xb8\x9d\x9b\x07\x0b\xeaL\x93\x82ppgaX<\x06A\x8eE\x8b&\x9c\x822y\x86\xcb\x90\x9a \x89\x01\x0c\x13g\xbe\x99a\n\x1c \x99}op\xed 09%\x03J\x87\x80\xc80F+\x83\x00\xc80\x92\x7fOS\x82\xdbl\x87\x99\\\x17|\xfe\x83\xde\xe7\xb1\xfa^\xab\x80E\xb1\xa4\xf3\xfc\x0f\x19\xdf\xa3\xf4\x9b\xe5\xbf\x8f\xd48\xbdL\x9e\xc5/\xe53h~\xc7\xdd)\xc6\xed\x99J\xf7\x8f\xfe]\xa70l\xaff\x13\x00\x9f\x0f^\xba\xde\x99\xefq\xaft\xbf\xa2\xae\xb9k\xcd\x84\xc1\x18!\xcd\x89\x83s\xbc\x06\xee\x91\x00\xf95\xe3\x909\x13Ji\x17\xd6\xdc\xab\xa3\xf0\xa8\x11\xb0\x99C2k?%\x80A\xde\xcb;\xcd\x1b\x99\x9b)\xb8\xd7\xea\xab\xbc\xb8C\xabUl\xb0\xd8_Y+!D\xc9o\x8c|\x9f)`.2j\xbd\x12\x19_\xd09-#un\x1e\x83\xdb!%\xe7H\xff'-\xe9>\x05H\xe7\xafd\xb8\xb1\x91w\x9e\xd5\x11\xea\xf3Q}\re\x10w\xf0jh\"/\x81\xebN\xe8\r\xb3\xcd\x0e\x04\x0b-:)\x05\xd3`\xab=\xd8\xa3\\ \x1a\xf1\xb6\xeei\\\xde\x9b\xb1\x13f\x81\x92GS\xd79\xae\xbd\x05K\xcf_\x83\x9b\xab\x0e\xb0\x7fr\\\xf5\xef\xf9\xbd(\xdd\xed\xe2)\xf6qT\xc3D\x9e\xfd0\xfcE*A\x0eDp\xea\xd5;7-U\xf65\xe5\xf6\xde\x8b\xb48\x82;G<\x9d\xd3:\xba\xca\xe1\x84\xf3c\xd3\xb0\xf4\x9d\xd8\x05\xe7\xec\xe1j\xd0d\x94>\xf3l-\xbc\x1a\xc4.z\x94,R\xbbA\x12\x9a\xa8v\xbb\x1b\x92\n\xab\xb8\xba\xe2t\xe3\x8b\xe4\xff\xd3\x1a\x14\xcc\xc2\x8dt\x16\xe3F\xc7>\xfd6\x03\x15\xcc\xe1x}\xadk\x97\x84\xc4\xc8\xa1\x82\x88\x81\x8d\xd3s\xba\xbc\xeej\xb5\x85\xd5}\xf8\xf8\xfe\xec\xad\x1b\x88\x11\x91\x03\xfb}\xfe\\\xa6_\xf8\xbb\x92)\xc2\x84\x84oP\xeb\xc8"

stat_files = [
    "stats/stats.pickle.xz",
    "stats/reads/reads.csv.bz2",
    "stats/reads/reads.png",
    "stats/HiCPairs/HiCPairs.csv.bz2",
    "stats/HiCPairs/HiCPairs.png",
    "stats/HiCPairs/intra_inter_sequence_valid_pairs.csv.bz2",
    "stats/HiCPairs/intra_inter_sequence_valid_pairs.png",
    "stats/sequences/per_sequence/per_sequence.csv.bz2",
    "stats/sequences/per_sequence/invalid_pairs.png",
    "stats/sequences/per_sequence/total_invalid_pairs.png",
    "stats/sequences/per_sequence/restriction_sites.png",
    "stats/sequences/per_sequence/inter_sequence_indexes.csv.bz2",
    "stats/sequences/per_sequence/inter_sequence_indexes.png",
    "stats/sequences/valid_pairs/valid_pairs.csv.bz2",
    "stats/sequences/valid_pairs/valid_pairs_raw_counts.png",
    "stats/sequences/valid_pairs/valid_pairs_norm_counts.png",
    "stats/sequences/valid_pairs/total_valid_pairs_raw_count.png",
    "stats/sequences/valid_pairs/total_valid_pairs_norm_counts.png",
    "stats/restriction_sites/restriction_sites.csv.bz2",
    "stats/restriction_sites/valid_pairs/pair_type.png",
    "stats/restriction_sites/valid_pairs/read_dir.png",
    "stats/restriction_sites/valid_pairs/fr_rf_bivariate.png",
    "stats/restriction_sites/invalid_pairs/pair_type.png",
    "stats/restriction_sites/invalid_pairs/read_dir.png",
    "stats/restriction_sites/invalid_pairs/sc_de_bivariate.png",
    "stats/intra_sequence_pair_separations/intra_sequence_pair_separations.csv.bz2",
    "stats/intra_sequence_pair_separations/histograms_plus_models/histograms_plus_models.csv.bz2",
    "stats/intra_sequence_pair_separations/histograms_plus_models/pair_separation.png",
    "stats/intra_sequence_pair_separations/histograms_plus_models/pair_separation_by_type.png",
    "stats/intra_sequence_pair_separations/histograms_plus_models/model_fits_by_pair_type.png",
    "stats/intra_sequence_pair_separations/histograms_plus_models/model_fits.png",
    "stats/intra_sequence_pair_separations/model_F_tests/model_F_tests.csv.bz2",
    "stats/intra_sequence_pair_separations/model_F_tests/F_tests_by_pair_type.png",
    "stats/intra_sequence_pair_separations/model_F_tests/F_tests.png",
    "stats/intra_sequence_fragment_separations/intra_sequence_fragment_separations.csv.bz2",
    "stats/intra_sequence_fragment_separations/histograms_plus_models/histograms_plus_models.csv.bz2",
    "stats/intra_sequence_fragment_separations/histograms_plus_models/pair_separation.png",
    "stats/intra_sequence_fragment_separations/histograms_plus_models/pair_separation_by_type.png",
    "stats/intra_sequence_fragment_separations/histograms_plus_models/model_fits_by_pair_type.png",
    "stats/intra_sequence_fragment_separations/histograms_plus_models/model_fits.png",
    "stats/intra_sequence_fragment_separations/model_F_tests/model_F_tests.csv.bz2",
    "stats/intra_sequence_fragment_separations/model_F_tests/F_tests_by_pair_type.png",
    "stats/intra_sequence_fragment_separations/model_F_tests/F_tests.png",
]


def test_hiline(tmpdir):
    os.chdir(tmpdir)

    with Popen("HiLine", stdout=PIPE, stderr=STDOUT) as process:
        process.communicate()

    with open("reads.cram", "wb") as file:
        file.write(bz2.decompress(reads))
    with open("ref.fa.gz", "wb") as file:
        file.write(bz2.decompress(ref))

    pipeline = Pipeline()
    pipeline.reference = "ref.fa.gz"
    pipeline.restriction_sites = "Bsp19I"
    pipeline.threads = 4
    pipeline.min_mapq = 10
    pipeline.register_bwa_alignment_sam_reads(reads="reads.cram", mark_dups=True)
    for handle in pipeline.output.all_reads:
        pipeline.register_output_file(file_name="all.cram", sort=True, handle=handle)
    pipeline.save_stats_path = "stats"
    pipeline.run()

    read, write = os.pipe()

    def thread_fn():
        with os.fdopen(write, "wb") as file:
            file.write(bz2.decompress(all_cram))

    thread = Thread(target=thread_fn)
    thread.start()

    with Popen(
        "samtools view --no-PG -h all.cram".split(), stdout=PIPE, stderr=STDOUT
    ) as process1, Popen(
        "samtools view --reference ref.fa.gz --no-PG -h -".split(),
        stdin=read,
        stdout=PIPE,
        stderr=STDOUT,
    ) as process2:
        replaceVN = True
        for i, (line1, line2) in enumerate(zip(process1.stdout, process2.stdout)):
            if 0 < i < 5:
                line1 = line1.decode("utf-8", errors="ignore").split("\t")
                line1 = ("\t".join(line1[:-1]) + "\n").encode("utf-8")
                line2 = line2.decode("utf-8", errors="ignore").split("\t")
                line2 = ("\t".join(line2[:-1]) + "\n").encode("utf-8")

            if replaceVN:
                if (
                    "@PG".encode("utf-8") in line2
                    and "ID:HiLine".encode("utf-8") in line2
                ):
                    replaceVN = False
                    line2 = re.sub(
                        r"VN:\d+\.\d+\.\d+",
                        "VN:{ver}".format(ver=version),
                        line2.decode("utf-8", errors="ignore"),
                    ).encode("utf-8")

            assert line1 == line2

    thread.join()

    test_stat_files = []
    for root, _, filenames in os.walk("stats"):
        for file in filenames:
            test_stat_files.append(os.path.join(root, file))
    assert len(test_stat_files) == len(stat_files)
    for file in test_stat_files:
        assert file in stat_files


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_hiline(tmpdir)
