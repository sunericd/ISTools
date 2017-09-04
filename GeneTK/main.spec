a = Analysis(['main.py'],
             pathex=['C:\\Users\\edsun\\Desktop\\IntegratedSciences\\ISTools\\GeneTK'],
             datas = [('data/restriction_sites2.csv', 'Data'), ('data/logo.gif', 'Data')],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)


pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='GeneTK.exe',
          debug=False,
          strip=None,
          upx=True,
          console=False , icon='C:\\Users\\edsun\\Desktop\\IntegratedSciences\\ISTools\\GeneTK\\data\\logo.ico')