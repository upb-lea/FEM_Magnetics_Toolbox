

arr_num_checks = [self.N95_checkBox.isChecked(), self.N97_checkBox.isChecked(), self.N87_checkBox.isChecked()]
# arr_choice_block = [self.mat_choice_1.setText(""), self.mat_choice_2.setText(), self.mat_choice_3.setText()]

for i in arr_num_checks:
    if arr_num_checks[i]:
        self.mat_choice_1.setText("N95")

    if arr_num_checks[1]:
        self.mat_choice_1.setText("N97")
    if arr_num_checks[2]:
        self.mat_choice_1.setText("N87")

    elif arr_num_checks[1]:
        self.mat_choice_1.setText("N97")
        if arr_num_checks[2]:
            self.mat_choice_2.setText("N87")
    elif arr_num_checks[2]:
        self.mat_choice_1.setText("N87")

    elif not arr_num_checks[0]:
        self.mat_choice_1.setText("")
    elif not arr_num_checks[1]:
        self.mat_choice_2.setText("")
    elif not arr_num_checks[2]:
        self.mat_choice_3.setText("")


