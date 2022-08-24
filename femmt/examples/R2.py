from PyQt4.QtGui import *
 model = QStandardItemModel()
 item = QStandardItem("Item")
 item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
 item.setData(QVariant(Qt.Checked), Qt.CheckStateRole)
 model.appendRow(item)
 view = QListView()
 view.setModel(model)
 view.show()