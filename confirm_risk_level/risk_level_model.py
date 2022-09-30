"""
@File ：fetch_gene_part.py
@Author ：wu_leanne
@Date ：2022/9/30
@Desc : This module is used for confirming risk level
@Email : wu_leanne@163.com
"""

import re


class ChinaLevel:
    def __init__(self, instestines_obj):
        self.instestines_obj = instestines_obj
        self.xianliu_count = 0
        self.level = "——"
        self.classify_level()

    def classify_level(self):
        """
        :return:
        """
        # 判断是否复核
        if self.need_check():
            return "需复核"
        # 判断是否为7
        if self.classify_level_7():
            self.level = 7
            return
        if self.classify_level_6():
            self.level = 6
            return
        if self.classify_level_5():
            self.level = 5
            return
        if self.classify_level_4():
            self.level = 4
            return
        if self.classify_level_3():
            self.level = 3
            return
        if self.classify_level_2():
            self.level = 2
            return
        if self.classify_level_1():
            self.level = 1
            return
        self.level = "——"
        return

    def need_check(self):
        return False

    def classify_level_7(self):
        count = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤|增生.{0,3}息肉.{0,3}伴.{0,6}腺瘤", item['pathology']) or re.search("腺瘤性息肉",
                                                                                         item['pathology']):
                count += item['number']
        self.xianliu_count = count
        if count > 10:
            return True
        return False

    def classify_level_6(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search("上皮内瘤变|异型增生", item['pathology']) and re.search("锯齿状", item['pathology']):
                return True
            if item['max_size'] >= 1 and re.search("锯齿状", item['pathology']) and (
                    not re.search("传统", item['pathology'])):
                return True
            if re.search('锯齿状腺瘤', item['pathology']):
                return True
            if re.search('传统.{0,6}锯齿状', item['pathology']):
                return True
        return False

    def classify_level_5(self):
        for item in self.instestines_obj.trans_cal_result:
            if item['max_size'] >= 1 and re.search('腺瘤', item['pathology']):
                return True
            if re.search("绒毛.{0,3}状", item['pathology']):
                return True
            if re.search("腺瘤", item['pathology']) and re.search("高级别", item['pathology']) and re.search("上皮内瘤|异型增生",
                                                                                                        item[
                                                                                                            'pathology']):
                return True
        return False

    def classify_level_4(self):
        if 3 <= self.xianliu_count < 10:
            return True
        return False

    def classify_level_3(self):
        if 1 <= self.xianliu_count <= 2:
            return True
        return False

    def classify_level_2(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿", item['pathology']) and (not re.search('上皮内瘤|异型增生', item['pathology'])) and item[
                'max_size'] < 1:
                return True
        return False

    def classify_level_1(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search('增生.{0,3}息肉', item['pathology']) and item['max_size'] >= 1:
                return False
        return True


class JapanLevel:
    def __init__(self, instestines_obj):
        self.instestines_obj = instestines_obj
        self.xianliu_count = 0
        self.level = "——"
        self.classify_level()

    def classify_level(self):
        """
        # 定级
        :return:
        """
        # 判断是否复核
        if self.need_check():
            return "需复核"
        # 判断是否为7
        if self.classify_level_5():
            self.level = 5
            return
        if self.classify_level_4():
            self.level = 4
            return
        if self.classify_level_3():
            self.level = 3
            return
        if self.classify_level_2():
            self.level = 2
            return
        self.level = 1
        return

    def need_check(self):
        return False

    def classify_level_5(self):
        count = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤", item['pathology']):
                count += item['number']
        self.xianliu_count = count
        if self.xianliu_count >= 10:
            return True
        for item in self.instestines_obj.trans_cal_result:
            if item['max_size'] >= 1 and re.search('腺瘤', item['pathology']):
                return True
            if re.search("绒毛.{0,3}状", item['pathology']):
                return True
            if re.search("高级别", item['pathology']) and re.search("上皮内瘤|异型增生", item[
                'pathology']):
                return True
        return False

    def classify_level_4(self):
        if 3 <= self.xianliu_count < 10:
            return True
        return False

    def classify_level_3(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿状", item['pathology']) and (re.search("无蒂", item['pathology'])):
                return True
        return False

    def classify_level_2(self):
        if 1 <= self.xianliu_count <= 2:
            return True
        return False


class EuropeLevel:
    def __init__(self, instestines_obj):
        self.instestines_obj = instestines_obj
        self.xianliu_count = 0
        self.level = "——"
        self.classify_level()

    def classify_level(self):
        """
        # 定级
        :return:
        """
        # 判断是否复核
        if self.need_check():
            return "需复核"
        # 判断是否为7
        if self.classify_level_5():
            self.level = 5
            return
        if self.classify_level_4():
            self.level = 4
            return
        if self.classify_level_3():
            self.level = 3
            return
        if self.classify_level_2():
            self.level = 2
            return
        # if self.classify_level_1():
        #     self.level = 1
        #     return
        self.level = 1
        return

    def need_check(self):
        return False

    def classify_level_5(self):
        count = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤|低级别.{0,6}上皮内瘤变|轻中度.{0,6}异型增生", item['pathology']):
                count += item['number']
        self.xianliu_count = count
        if self.xianliu_count >= 10:
            return True
        return False

    def classify_level_4(self):
        for item in self.instestines_obj.trans_cal_result:
            if item['max_size'] >= 2:
                return True
        return False

    def classify_level_3(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿状", item['pathology']) and re.search("传统", item['pathology']):
                return True
        for item in self.instestines_obj.trans_cal_result:
            if item['max_size'] >= 1 and re.search('锯齿状', item['pathology']):
                return True
            if re.search('锯齿状', item['pathology']) and (re.search('高级别|腺瘤|上皮内瘤变|异型增生', item['pathology'])):
                return True
        return False

    def classify_level_2(self):
        if self.xianliu_count >= 5:
            return True
        for item in self.instestines_obj.trans_cal_result:
            print(item['max_size'])
            if re.search("腺瘤|低级别.{0,6}上皮内瘤变|轻中度.{0,6}异型增生", item['pathology']) and item['max_size'] >= 1:
                return True
            if re.search("腺瘤|低级别.{0,6}上皮内瘤变|轻中度.{0,6}异型增生", item['pathology']) and re.search("高度|高级别",
                                                                                             item['pathology']):
                return True
        return False

    def classify_level_1(self):
        if 1 <= self.xianliu_count < 5:
            return True
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤|低级别.{0,6}上皮内瘤变|轻中度.{0,6}异型增生", item['pathology']) and item['max_size'] < 1:
                return True
            if re.search("锯齿状", item['pathology']) and item['max_size'] < 1 and (
                    not re.search('高级别|腺瘤|上皮内瘤变|异型增生', item['pathology'])):
                return True
        return False


class BritainLevel:
    def __init__(self, instestines_obj):
        self.instestines_obj = instestines_obj
        self.xianliu_count = 0
        self.level = "——"
        self.classify_level()

    def classify_level(self):
        """
        # 定级
        :return:
        """
        # 判断是否复核
        if self.need_check():
            return "需复核"
        if self.classify_level_4():
            self.level = 4
            return
        if self.classify_level_3():
            self.level = 3
            return
        if self.classify_level_2():
            self.level = 2
            return
        self.level = 1
        return

    def need_check(self):
        return False

    def classify_level_4(self):
        for item in self.instestines_obj.trans_cal_result:
            if item['max_size'] >= 2 and re.search("无蒂", item['pathology']):
                return True
        return False

    def classify_level_3(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿状", item['pathology']) or re.search("腺瘤|异型增生", item['pathology']):
                self.xianliu_count += item['number']
        if self.xianliu_count >= 5:
            return True
        for item in self.instestines_obj.trans_cal_result:
            if re.search('无蒂', item['pathology']):
                return True
        return False

    def classify_level_2(self):
        if self.xianliu_count >= 2:
            for item in self.instestines_obj.trans_cal_result:
                if re.search("锯齿状|腺瘤", item['pathology']) and item['max_size'] >= 1:
                    print(1)
                    return True
                if re.search("锯齿状", item['pathology']) and re.search("异型增生", item['pathology']):
                    print(2)
                    return True
                if re.search("腺瘤", item['pathology']) and re.search("锯齿状", item['pathology']):
                    print(3)
                    return True
                if re.search("高级别", item['pathology']):
                    print(4)
                    return True
        return False


class AmericaLevel:
    def __init__(self, instestines_obj):
        self.instestines_obj = instestines_obj
        self.xianliu_count = 0
        self.level = "——"
        self.classify_level()

    def classify_level(self):
        """
        # 定级
        :return:
        """
        # 判断是否复核
        if self.need_check():
            return "需复核"
        if self.classify_level_7():
            self.level = 7
            return
        if self.classify_level_6():
            self.level = 6
            return
        if self.classify_level_5():
            self.level = 5
            return
        if self.classify_level_4():
            self.level = 4
            return
        if self.classify_level_3():
            self.level = 3
            return
        if self.classify_level_2():
            self.level = 2
            return
        self.level = 1
        return

    def need_check(self):
        return False

    def classify_level_7(self):
        for item in self.instestines_obj.trans_cal_result:
            if item['max_size'] >= 2 and re.search("锯齿状", item['pathology']):
                return True
        return False

    def classify_level_6(self):
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤", item['pathology']):
                self.xianliu_count += item['number']
        if self.xianliu_count >= 10:
            return True
        return False

    def classify_level_5(self):
        xianliu_count1 = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤", item['pathology']) and item['max_size'] < 1:
                xianliu_count1 += item['number']
        if 5 <= xianliu_count1 <= 10:
            return True
        juci_count = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿状", item['pathology']):
                juci_count += item['number']
        if 5 <= juci_count <= 10:
            return True
        for item in self.instestines_obj.trans_cal_result:
            if re.search("腺瘤|锯齿状", item['pathology']) and item['max_size'] >= 1:
                return True
            if re.search("高级别|重度", item['pathology']) and re.search("上皮内瘤变|异型增生|腺瘤", item['pathology']):
                return True
            if re.search("绒毛.{0,3}状", item['pathology']) and re.search("腺瘤|异型增生", item['pathology']):
                return True
            if re.search("传统", item['pathology']) and re.search("腺瘤|异型增生", item['pathology']) and re.search("锯齿状", item[
                'pathology']):
                return True
            if re.search("锯齿状", item['pathology']) and re.search("腺瘤|异型增生", item['pathology']):
                return True
        return False

    def classify_level_4(self):
        if 3 <= self.xianliu_count <= 4:
            return True
        juci_count = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿状", item['pathology']):
                juci_count += item['number']
        if 3 <= juci_count <= 4:
            return True
        for item in self.instestines_obj.trans_cal_result:
            if re.search("增生性.{0,6}息肉", item['pathology']) and item['max_size'] >= 1:
                return True
        return False

    def classify_level_3(self):
        juci_count = 0
        for item in self.instestines_obj.trans_cal_result:
            if re.search("锯齿状", item['pathology']):
                juci_count += item['number']
        if 1 <= juci_count <= 2:
            return True
        return False

    def classify_level_2(self):
        if 1 <= self.xianliu_count <= 2:
            return True
        return False
