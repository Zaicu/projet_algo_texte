import sys
import os
import time
import platform
import inspect
import threading
from PyQt5.QtWidgets import QApplication, QWidget, QTreeView, QFileSystemModel, QVBoxLayout, QGridLayout, QRadioButton, QPushButton, QLabel
from PyQt5.QtCore import QModelIndex
src_file_path = inspect.getfile(lambda: None)
root_dir = os.path.dirname(os.path.realpath(src_file_path))
dir = root_dir+r'\Test'
print(root_dir)
sys.path.append(root_dir)
from parsing import *

if platform.system() == "Windows": SEP = "\\"
else: SEP = "/"


class FileSystemView(QWidget):
    def __init__(self, dir_path):
        self.path = "blank"
        super().__init__()
        appWidth = 800
        appHeight = 300
        self.setWindowTitle('File System Viewer')
        self.setGeometry(300, 300, appWidth, appHeight)

        self.model = QFileSystemModel()
        self.model.setRootPath(dir_path)
        self.tree =  QTreeView()
        self.tree.setModel(self.model)
        self.tree.setRootIndex(self.model.index(dir))
        self.tree.setColumnWidth(0, 250)
        self.tree.setAlternatingRowColors(True)
        self.tree.setStyleSheet("background-color:slategray; alternate-background-color:lightslategray; border-radius:5px")
        layout = QVBoxLayout()
        layout.addWidget(self.tree)
        self.setLayout(layout)

        self.tree.clicked.connect(self.test)

    def test(self, signal):
        self.content=self.tree.model().fileName(signal)
        dir2= dir.replace('\\','/')+'/'
        path1 = self.tree.model().filePath(signal)
        path2 = path1.replace(dir2,'')
        if platform.system()=='Windows':
            self.path = path2.replace("/","\\")
        else :
            self.path = path2


class Menu(QWidget):
    def __init__(self,list):
        super().__init__()
        layout = QGridLayout()
        self.setLayout(layout)
        for i in range(len(list)):
            content = list.pop()
            radiobutton=QRadioButton(content)
            if i==0 :
                radiobutton.setChecked(True)
                self.content = content
            radiobutton.content = content
            radiobutton.toggled.connect(self.onClicked)
            layout.addWidget(radiobutton, i, 0)


    def onClicked(self):
        radioButton = self.sender()
        self.content = radioButton.content



class Button(QWidget):
    def __init__(self,menu1,menu2):
        super().__init__()


        entree = QPushButton('Parse', self)
        entree.setCheckable(True)

        entree.clicked.connect(self.parse)
        entree.setStyleSheet("margin: 1px; padding: 10px; \
                           background-color: lightgray; \
                           border-style: solid; \
                           border-radius: 4px; border-width: 3px; \
                           border-color: slategray;")
        layout = QVBoxLayout()
        layout.addWidget(entree)
        self.setLayout(layout)
    def parse(self) :
        if(tree.path != "blank"):
            logs.write(tree.path+ " - ",menu_regions.content)
        th = threading.Thread(target=init,args=([menu_regions.content],root_dir))
        th.start()
        #init(menu_regions.content,root_dir)


class Button_init(QWidget):
    def __init__(self,menu1,menu2):
        super().__init__()


        self.entree = QPushButton('Init', self)
        self.entree.setCheckable(True)

        self.entree.clicked.connect(self.initialisation)
        self.entree.setStyleSheet("margin: 1px; padding: 10px; \
                           background-color: lightgray; \
                           border-style: solid; \
                           border-radius: 4px; border-width: 3px; \
                           border-color: slategray;")
        layout = QVBoxLayout()
        layout.addWidget(self.entree)
        self.setLayout(layout)
    def initialisation(self) :

        logs.write('',"Initialisation")
        th_init.start()
        #init(dir=root_dir)


class Logs(QWidget):
    def __init__(self):
        super().__init__()
        #self.setGeometry(200, 250, 100, 50)
        self.text= QLabel()
        self.text.setStyleSheet("background-color:slategray; border-radius:5px")
        self.text.setText("Logs")
        self.text.list = ["Logs"]
        self.text.setWordWrap(True)
        layout = QVBoxLayout()
        layout.addWidget(self.text)
        self.setLayout(layout)
    def write(self, kingdom, region):
        if len(self.text.list) > 10 :
            del(self.text.list[1])
        add = "Parsing en cours : "+ kingdom + region
        self.text.list.append("\n" + add)
        self.text.setText(" ".join(self.text.list))


if __name__ == '__main__':
    th_init = threading.Thread(target=init,args=([],root_dir))
    if not os.path.isdir(root_dir + SEP+ "GENOME_REPORTS"):
        os.makedirs(root_dir + SEP+ "GENOME_REPORTS")
    if not os.path.isdir(root_dir + SEP + "Results"):
        os.makedirs(root_dir + SEP + "Results")

    root_dir = os.path.dirname(os.path.realpath(src_file_path))
    dir = root_dir+SEP + 'Results'
    sys.path.append(root_dir)
    app = QApplication([])
    grid = QGridLayout()
    win = QWidget()
    tree = FileSystemView(dir)
    #menu_kingdom = Menu(["Eucaryota", "Bacteria", "Archaea","Virus","Plasmides", "Organelle"])
    menu_regions = Menu(["5'UTR","3'UTR","tRNA","telomere","rRNA","ncRNA","mobile_element","intron","centromere","CDS"])
    logs = Logs()
    button = Button(tree, menu_regions)
    button_init = Button_init(tree, menu_regions)

    grid.addWidget(tree,1,1,5,1)
    grid.addWidget(logs,3,2,2,2)
    #grid.addWidget(menu_kingdom,1,2)
    grid.addWidget(menu_regions,1,3)
    grid.addWidget(button_init,2,2)
    grid.addWidget(button,2,3)


    win.setLayout(grid)
    win.setWindowTitle("Logiciel Bioinformatique")
    win.show()
    sys.exit(app.exec_())
