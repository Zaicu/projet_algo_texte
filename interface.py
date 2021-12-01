import sys
import os
import time
import platform
import inspect
from PyQt5.QtWidgets import QApplication, QWidget, QTreeView, QFileSystemModel, QVBoxLayout, QGridLayout, QRadioButton, QPushButton, QLabel, QDialog, QProgressDialog, QButtonGroup, QLineEdit
from PyQt5.QtCore import QModelIndex, QRunnable, QThreadPool, QThread, pyqtSignal, Qt, QThread, QObject

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
        appWidth = 1000
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
        self.regions_group = QButtonGroup(self)

        self.setLayout(layout)
        self.nb = len(list)
        for i in range(len(list)):
            content = list.pop()
            radiobutton=QRadioButton(content)
            if i==0 :
                radiobutton.setChecked(True)
                self.content = [content]
            radiobutton.content = content
            radiobutton.toggled.connect(self.onClicked)
            layout.addWidget(radiobutton, i, 0)
            self.regions_group.addButton(radiobutton)
        self.regions_group.setExclusive(False)


    def onClicked(self):
        radioButton = self.sender()
        if self.sender().isChecked():
            self.content+=[radioButton.content]
        else :
            self.content.remove(radioButton.content)

    def select_all(self):
        for button in self.regions_group.buttons():
            button.setChecked(True)

    def unselect_all(self):
        for button in self.regions_group.buttons():
            button.setChecked(False)


class Button_regions(QWidget):
    def __init__(self,menu):
        super().__init__()

        self.menu = menu
        self.total = menu.nb
        entree = QPushButton('Select/Unselect all regions', self)
        entree.setCheckable(True)

        entree.clicked.connect(self.select)
        #entree.setStyleSheet("margin: 1px; padding: 10px; \
        #                   background-color: lightgray; \
         #                  border-style: solid; \
          #                 border-radius: 4px; border-width: 3px; \
           #                border-color: slategray;")
        layout = QVBoxLayout()
        entree.setMaximumSize(200,30)
        layout.addWidget(entree)
        self.setLayout(layout)


    def select(self) :
        self.selected=self.menu.content
        if len(self.selected)==self.total:
            self.menu.unselect_all()
        else :
            self.menu.select_all()



class Button_init(QWidget):

    def __init__(self,menu1,menu2,text):
        super().__init__()


        self.entree = QPushButton('Parse', self)
        self.entree.setCheckable(True)
        self.text = text.qle

        self.entree.clicked.connect(self.parse)
        self.entree.setStyleSheet("margin: 1px; padding: 10px; \
                           background-color: lightgray; \
                           border-style: solid; \
                           border-radius: 4px; border-width: 3px; \
                           border-color: slategray;")
        layout = QVBoxLayout()
        layout.addWidget(self.entree)
        self.setLayout(layout)
        self.threadpool = QThreadPool()

    def process_results(self,tuple):
        self.list_init=tuple
        logs.write("You can now parse")
        #print(self.list_init)

    def initialisation(self):
        logs.write("Please wait during the download")
        init_worker = Worker()
        self.threadpool.start(init_worker)
        #self.threadpool.waitForDone()
        init_worker.signals.result.connect(self.process_results)

    def parse(self) :
        if(tree.path != "blank"):
            to_write =''
            for i in range(len(menu_regions.content)):
                if i==len(menu_regions.content)-1:
                    to_write+=menu_regions.content[i]
                else:
                    to_write+=menu_regions.content[i]+", "
            if self.text.text():
                if len(menu_regions.content)==0:
                    to_write+=self.text.text()
                else:
                    to_write+=", "+self.text.text()
            logs.write_parse(tree.path+ " - ",to_write)
            parse_worker = Worker_parse()
            if self.text.text():
                parse_worker.set_param(self.list_init,tree.content,menu_regions.content+[self.text.text()])
            else:
                parse_worker.set_param(self.list_init,tree.content,menu_regions.content)
            self.threadpool.start(parse_worker)
            parse_worker.signal.finished.connect(self.finished)

        else:
            to_write =''
            for i in range(len(menu_regions.content)):
                if i==len(menu_regions.content)-1:
                    to_write+=menu_regions.content[i]
                else:
                    to_write+=menu_regions.content[i]+", "
            if self.text.text():
                if len(menu_regions.content)==0:
                    to_write+=self.text.text()
                else:
                    to_write+=", "+self.text.text()
            logs.write_parse("Archaea, Bacteria, Eukaryota, Viruses"+ " - ",to_write)
            parse_worker = Worker_parse()
            if self.text.text():
                parse_worker.set_param(self.list_init,os.path.join("Results"),menu_regions.content+[self.text.text()])
            else:
                parse_worker.set_param(self.list_init,os.path.join("Results"),menu_regions.content)

            self.threadpool.start(parse_worker)
            parse_worker.signal.finished.connect(self.finished)
            #self.threadpool.waitForDone()
            #associate(self.list_init[0], self.list_init[1], self.list_init[2],tree.content,menu_regions.content)
            #associate([menu_regions.content,tree.content])


    def finished(self):
        logs.write("Finished parsing")

class WorkerSignals(QObject):
    result = pyqtSignal(list)

class Worker(QRunnable):
    def __init__(self):
        super(Worker, self).__init__()
        self.signals = WorkerSignals()

    def run(self):
        self.ids, self.path, self.dates=init(logs, prgss)
        self.list= [self.ids, self.path, self.dates]
        time.sleep(5)
        self.signals.result.emit(self.list)


class Signals(QObject):
    finished = pyqtSignal()

class Worker_parse(QRunnable):
    def __init__(self):
        super(Worker_parse, self).__init__()
        self.signal = Signals()

    def set_param(self,list_init, tree, regions):
        self.list_init=list_init
        self.tree=tree
        self.regions=regions

    def run(self):
        associate(self.list_init[0], self.list_init[1], self.list_init[2],self.tree,self.regions, prgss_parsing)
        time.sleep(5)
        self.signal.finished.emit()

class FreeText(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):

        hbox = QVBoxLayout(self)
        self.qle = QLineEdit(self)

        self.label = QLabel("Free region text :")
        hbox.addWidget(self.label)

        hbox.addSpacing(20)
        hbox.addWidget(self.qle)
        self.setLayout(hbox)

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
    def write_parse(self, kingdom, region):
        if len(self.text.list) > 10 :
            del(self.text.list[1])
        add = "Ongoing parsing : "+ kingdom + region
        self.text.list.append("\n" + add)
        self.text.setText(" ".join(self.text.list))
    def write(self, text):
        if len(self.text.list) > 10 :
            del(self.text.list[1])
        self.text.list.append("\n" + text)
        self.text.setText(" ".join(self.text.list))
    def prg(self,total):
        print("creation ",total)
        #prgss.create(total)
        #self.prgss=ProgressBar(total)
        #print("created")
    def update(self, value,total):
        prgss.signal_update.emit(int(value/total)*100)
        print(value,total)
        if(value==total):
            print("value=total")
            #wait(1)
            prgss.close()
            prgss.setVisible(False)



class ProgressBar(QProgressDialog):
    signal_update=pyqtSignal(int)

    def __init__(self, max):
        #print("init")
        super().__init__()
        self.setMinimumDuration(0) # Sets how long the loop should last before progress bar is shown (in miliseconds)
        self.setWindowTitle("Creating directories")
        self.setModal(True)
        #self.setWindowModality(Qt::WindowModal)

        self.setValue(0)
        self.setMinimum(0)
        self.setMaximum(max)
        self.signal_update.connect(self.update)
        #self.setAutoReset(True)
        #self.setAutoClose(True)
        #self.setAttribute(Qt::WA_DeleteOnClose, true);
        self.resize(300,100)
        self.show()

    def update(self, value):
        self.setValue(value)
        self.show()
        if value==100:
            self.setVisible(False)
            self.close()

    #def close(self):
    #    print("allo",self.wasCanceled())
    #    self.close()
    #    self.setVisible(False)

class ProgressBar_parsing(QProgressDialog):
    signal_update=pyqtSignal(int)

    def __init__(self, max):
        super().__init__()
        self.setMinimumDuration(0) # Sets how long the loop should last before progress bar is shown (in miliseconds)
        self.setWindowTitle("Parsing")
        self.setModal(True)

        self.setValue(0)
        self.setMinimum(0)
        self.setMaximum(max)
        self.signal_update.connect(self.update)
        self.resize(300,100)
        self.show()

    def update(self, value):
        self.setValue(value)
        self.show()
        if value==100:
            self.setVisible(False)
            #self.close()

if __name__ == '__main__':
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
    prgss = ProgressBar(100)
    prgss.setVisible(False)

    prgss_parsing = ProgressBar(100)
    prgss_parsing.setVisible(False)

    freetext = FreeText()

    button_select = Button_regions(menu_regions)
    button_init = Button_init(tree, menu_regions, freetext)
    label_menu = QLabel("Select one or several region(s) :")

    grid.addWidget(tree,1,1,7,1)
    grid.addWidget(logs,6,2,2,2)
    #grid.addWidget(menu_kingdom,1,2)
    grid.addWidget(menu_regions,3,2,1,2)
    grid.addWidget(button_init,5,2,1,2)
    grid.addWidget(button_select,2,2,1,1)
    grid.addWidget(freetext,4,2,1,1)
    grid.addWidget(label_menu,1,2,1,1)


    win.setLayout(grid)
    win.setGeometry(300,300,1000,300)
    win.setWindowTitle("Logiciel Bioinformatique")
    win.show()
    button_init.initialisation()

    sys.exit(app.exec_())
