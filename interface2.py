import sys
import os
import platform
import inspect
from PyQt5.QtWidgets import QApplication, QWidget, QTreeView, QFileSystemModel, QVBoxLayout, QGridLayout, QRadioButton, QPushButton, QLabel
from PyQt5.QtCore import QModelIndex
from parsing2 import *
src_file_path = inspect.getfile(lambda: None)


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
		self.tree.setRootIndex(self.model.index(dirPath))
		self.tree.setColumnWidth(0, 250)
		self.tree.setAlternatingRowColors(True)

		layout = QVBoxLayout()
		layout.addWidget(self.tree)
		self.setLayout(layout)

		self.tree.clicked.connect(self.test)

    def test(self, signal):
        self.content=self.tree.model().fileName(signal)
        dirPath2= dirPath.replace('\\','/')
        path1 = self.tree.model().filePath(signal)
        path2 = path1.replace(dirPath2,'')
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
		layout = QVBoxLayout()
		layout.addWidget(entree)
		self.setLayout(layout)
	def parse(self) :
	    if(tree.path != "blank"):
            logs.write(tree.path,menu_regions.content)

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
		add = "Parsing en cours : "+ kingdom + " - "+ region
		self.text.list.append("\n" + add)
		self.text.setText(" ".join(self.text.list))


if __name__ == '__main__':
	app = QApplication([])
	grid = QGridLayout()
	win = QWidget()

	str = os.path.dirname(os.path.realpath(src_file_path))
	dirPath = str+r'\Test'
    sys.path.append(str)
	tree = FileSystemView(dirPath)
	#menu_kingdom = Menu(["Eucaryota", "Bacteria", "Archaea","Virus","Plasmides", "Organelle"])
	menu_regions = Menu(["CDS", "Intron", "Télomère"])
	logs = Logs()
	button = Button(tree, menu_regions)

	grid.addWidget(tree,1,1,5,1)
	grid.addWidget(logs,3,2,2,2)
	#grid.addWidget(menu_kingdom,1,2)
	grid.addWidget(menu_regions,1,3)
	grid.addWidget(button,2,3)


	win.setLayout(grid)
	win.setWindowTitle("Logiciel Bioinformatique")
	win.show()
	sys.exit(app.exec_())
