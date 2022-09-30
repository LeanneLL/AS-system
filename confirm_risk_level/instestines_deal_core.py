"""
@File ：fetch_gene_part.py
@Author ：wu_leanne
@Date ：2022/9/30
@Desc : This module is the confirm_risk_level of processing data
@Email : wu_leanne@163.com
"""
import re
from confirm_risk_level.ai_model import AiSite, AiNumber, AiSize, AiBiopsy
from confirm_risk_level.pathlogy_model import PathlogyRegDesc
from confirm_risk_level.risk_level_model import ChinaLevel, JapanLevel, EuropeLevel, BritainLevel, AmericaLevel


class IntestinesObject:
    """
    confirm_risk_level of processing  model
    """

    def __init__(self, check_id="", microscopic_desc="", microscopic_diagnosis="", pathology="", level_standard="China",
                 entity_list=None):
        self.level = 0
        self.ruzu = None
        self.check_id = check_id
        self.microscopic_desc = microscopic_desc if str(microscopic_desc) != "nan" else ""
        self.microscopic_diagnosis = microscopic_diagnosis if str(microscopic_diagnosis) != "nan" else ""
        self.pathology = pathology.replace("_x005f", "") if str(pathology) != "nan" else ""
        self.pathology = pathology.replace("_x005f", "")
        self.level_standard = level_standard
        # 部位类列表
        self.ai_site_list = entity_list[0]['polyp site'] if 'polyp site' in entity_list[0] else []
        # 数量类列表
        self.ai_number_list = entity_list[0]['polyp number'] if 'polyp number' in entity_list[0] else []
        # 大小类列表
        self.ai_size_list = entity_list[0]['polyp size'] if 'polyp size' in entity_list[0] else []
        # 送检类列表
        self.ai_songjian_list = entity_list[0]['polyp pathology'] if 'polyp pathology' in entity_list[0] else []
        # 部位对应距肛门距离范围
        self.site_len_dict = {
            "直肠": [0, 17],
            "乙状结肠": [15, 40],
            "降结肠": [30, 50],
            "横结肠": [45, 65],
            "升结肠": [55, 1000],
        }
        # 部位特殊写法
        self.special_site_dict = {
            "回盲": "回盲部|回盲瓣|阑尾开口|盲肠",
            "直肠": "近肛门口|紧邻肛缘"
        }
        self.ai_obj_list = []  # 串起来
        self.classifed_ai_keys = []  # 串起来分类后
        self.reg_site_ai_keys_dict = {}  # 与病理匹配后的
        self.reg_site_ai_obj_list_dict = {}  # 过滤送检后的
        self.trans_cal_result = []  # 转化后计算
        self.xianliu_count = 0  # 管状腺瘤数量
        # 生成病理类
        self.pathology_obj = PathlogyRegDesc(self.pathology)
        # 处理数据
        self.deal_data()

    def deal_data(self):
        # 预处理,将不规范的去除
        self.pre_deal()
        # 将所有关键字串起来,如果第一个不是部位则添加部位
        self.join_ai_keys()
        # 预测可能的部位
        self.predict_possible_site()
        # 再次将所有关键字串起来
        self.join_ai_keys()
        # 将串起来的部位分类
        self.classify_ai_keys()
        # 将部位和病理匹配
        self.match_site_with_pathology()
        # 清除有送病检但没写的
        self.clear_no_sonjian()
        # 转化格式并计算
        self.trans_cal()
        # 定级
        self.classify_level()
        # 是否入组
        self.classify_ruzu()

    def pre_deal(self):
        """
        # 预处理,将不规范的去除，如果传入字段不含有['text', 'start', 'end']三个键，则排除
        :return:
        """
        ai_site_list = []

        colums = ['text', 'start', 'end']
        for item in self.ai_site_list:
            valid = True
            for c in colums:
                if c not in item:
                    valid = False
                    break
            if valid:
                r_ = re.search("(\d+)cm以下", item['text'], flags=re.I)
                if r_:
                    item['text'] = re.sub("(\d+)", str(int(r_.group(1)) - 1), item['text'])
                # xx近回盲部则不算回盲部
                item['text'] = re.sub("近回盲部", "", item['text'])
                ai_site_list.append(item)
        self.ai_site_list = ai_site_list

        ai_number_list = []
        for item in self.ai_number_list:
            valid = True
            for c in colums:
                if c not in item:
                    valid = False
                    break
            if valid:
                ai_number_list.append(item)
        self.ai_number_list = ai_number_list

        ai_size_list = []
        for item in self.ai_size_list:
            if not re.search("\d", item['text']):
                continue
            valid = True
            for c in colums:
                if c not in item:
                    valid = False
                    break
            if valid:
                ai_size_list.append(item)
        self.ai_size_list = ai_size_list

        ai_songjian_list = []
        for item in self.ai_songjian_list:
            valid = True
            for c in colums:
                if c not in item:
                    valid = False
                    break
            if valid:
                ai_songjian_list.append(item)
        self.ai_songjian_list = ai_songjian_list
        # 将所有信息转化为对象
        self.ai_site_list = [AiSite(ai_keys['text'], ai_keys['start'], ai_keys['end']) for ai_keys in self.ai_site_list]
        self.ai_number_list = [AiNumber(ai_keys['text'], ai_keys['start'], ai_keys['end']) for ai_keys in
                               self.ai_number_list]
        self.ai_size_list = [AiSize(ai_keys['text'], ai_keys['start'], ai_keys['end']) for ai_keys in self.ai_size_list]
        self.ai_songjian_list = [AiBiopsy(ai_keys['text'], ai_keys['start'], ai_keys['end']) for ai_keys in
                                 self.ai_songjian_list]

    def join_ai_keys(self):
        """
        # 将所有关键字串起来，并排序
        :return:
        """
        ai_obj_list = []
        ai_obj_list += self.ai_site_list + self.ai_size_list + self.ai_number_list + self.ai_songjian_list
        self.ai_obj_list = sorted(ai_obj_list, key=lambda x: x.start)
        # 如果第一个对象不是部位对象，则人工添加一个对象
        if len(self.ai_obj_list) > 0 and self.ai_obj_list[0].type != "site":
            tmp_ai_site = AiSite("人工添加部位", 0, 1)
            self.ai_site_list.append(tmp_ai_site)

    def predict_possible_site(self):
        """
        # 预测可能的部位
        :return:
        """
        site_list = [
            "回盲部|回盲瓣|阑尾开口|盲肠|回盲",
            "升结肠",
            "横结肠",
            "降结肠",
            "乙状结肠",
            "直肠|近肛门口|紧邻肛缘",
        ]
        self.ai_site_list = sorted(self.ai_site_list, key=lambda x: x.start)
        # 双指针方式
        for i in range(len(self.ai_site_list)):
            # 找前面的
            l = i
            possible_start = 0
            while l >= 0:
                find = False
                for j in range(len(site_list)):
                    site_reg = site_list[j]
                    if re.search(site_reg, self.ai_site_list[l].text):
                        possible_start = j
                        find = True
                        break
                if find:
                    break
                l -= 1
            # # 找后面的
            r = i
            possible_end = len(site_list)
            while r <= len(self.ai_site_list) - 1:
                find = False
                for j in range(len(site_list)):
                    site_reg = site_list[j]
                    if re.search(site_reg, self.ai_site_list[r].text):
                        possible_end = j
                        find = True
                        break
                if find:
                    break
                r += 1
            # 设置可能部位
            self.ai_site_list[i].possible_site = list(
                site_list[possible_start:possible_end]) if possible_start != possible_end else list(
                site_list[possible_start:possible_end + 1])
            if re.search("乙状结肠", self.ai_site_list[i].text):
                self.ai_site_list[i].possible_site.append("直肠")

    def classify_ai_keys(self):
        """
        # 将串起来的部位分类
        :return:
        """
        self.classifed_ai_keys = []
        # 滑动窗口算法
        l = 0
        r = 0
        not_site_count = 0  # 窗口中非部位关键字数量
        while r < len(self.ai_obj_list):
            if self.ai_obj_list[r].type != "site":
                not_site_count += 1
            else:
                if not_site_count != 0:
                    tmp = list(self.ai_obj_list[l:r])
                    self.classifed_ai_keys.append(tmp)
                    l = r
                    not_site_count = 0
            r += 1

        tmp = list(self.ai_obj_list[l:r])
        self.classifed_ai_keys.append(tmp)

    def match_site_with_pathology(self):
        """
        # 将部位和病理匹配
        :return:
        """
        self.reg_site_ai_keys_dict = {}  # 结果
        reg_site_list = self.pathology_obj.site_reg_list
        if "降结肠|脾曲|脾区" in reg_site_list:
            reg_site_list.insert(0, "降结肠|脾曲|脾区")
            reg_site_list = list(set(reg_site_list))
        if "横结肠|脾曲" in reg_site_list:
            reg_site_list.insert(0, "横结肠|脾曲")
            reg_site_list = list(set(reg_site_list))
        # 将 结肠 部位放到最后面进行匹配，避免前面就讲其他部位匹配进去了
        jiechang_index = []
        for i in range(len(reg_site_list)):
            if reg_site_list[i] == "结肠":
                jiechang_index.append(i)
        if jiechang_index:
            jiechang = reg_site_list[jiechang_index[0]]
            last = reg_site_list[-1]
            reg_site_list[-1] = jiechang
            reg_site_list[jiechang_index[0]] = last
        matched_ai_index = []  # 已经匹配过的
        # 第一次匹配，直接匹，按照部位和病理直接匹配
        for reg_site in reg_site_list:
            self.reg_site_ai_keys_dict[reg_site] = []
            if not re.search(".+(肝曲|脾曲)", reg_site):
                # 如果部位的形式不是带肝曲脾曲的，则正常匹配
                for i in range(len(self.classifed_ai_keys)):
                    if i in matched_ai_index:
                        continue
                    ai_obj_list = self.classifed_ai_keys[i]
                    ai_site_list = [ai_obj.text for ai_obj in ai_obj_list if ai_obj.type == "site"]
                    ai_site_strs = ",".join(ai_site_list)
                    if re.search(reg_site, ai_site_strs):
                        self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                        matched_ai_index.append(i)
            else:
                # 部位对象中是否含有肝曲|脾曲
                find_piqu_anqu = False
                reg_site_ = reg_site.split("|")[1]  # 脾曲/肝曲
                for i in range(len(self.classifed_ai_keys)):
                    if i in matched_ai_index:
                        continue
                    ai_obj_list = self.classifed_ai_keys[i]
                    ai_site_list = [ai_obj.text for ai_obj in ai_obj_list if ai_obj.type == "site"]
                    ai_site_strs = ",".join(ai_site_list)
                    if re.search(reg_site_, ai_site_strs):
                        # 有肝曲和脾曲直接匹配
                        self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                        matched_ai_index.append(i)
                        find_piqu_anqu = True
                # 如果没有肝曲和脾曲，则按部位匹配
                if not find_piqu_anqu:
                    for i in range(len(self.classifed_ai_keys)):
                        if i in matched_ai_index:
                            continue
                        ai_obj_list = self.classifed_ai_keys[i]
                        ai_site_list = [ai_obj.text for ai_obj in ai_obj_list if ai_obj.type == "site"]
                        ai_site_strs = ",".join(ai_site_list)
                        if re.search(reg_site, ai_site_strs):
                            self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                            matched_ai_index.append(i)

        # 第二次匹配, 匹配特殊写法
        for reg_site in reg_site_list:
            if reg_site not in self.special_site_dict:
                continue
            for j in range(len(self.classifed_ai_keys)):
                if j in matched_ai_index:
                    continue
                ai_obj_list = self.classifed_ai_keys[j]
                ai_site_list = [ai_obj.text for ai_obj in ai_obj_list if ai_obj.type == "site"]
                ai_site_strs = ",".join(ai_site_list)
                reg_site_ = self.special_site_dict[reg_site]
                if re.search(reg_site_, ai_site_strs, flags=re.I):
                    self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                    matched_ai_index.append(j)

        # 第三次匹配 距离匹配 病理中是部位，镜下所见中是距离
        for reg_site in reg_site_list:
            if reg_site not in self.site_len_dict:
                continue
            for j in range(len(self.classifed_ai_keys)):
                if j in matched_ai_index:
                    continue
                ai_obj_list = self.classifed_ai_keys[j]
                ai_site_list = [ai_obj.text for ai_obj in ai_obj_list if ai_obj.type == "site"]
                # 找到提前预测好的可能部位
                ai_possible_site_list = []
                for ai_obj in ai_obj_list:
                    if ai_obj.type == "site":
                        ai_possible_site_list += ai_obj.possible_site
                ai_possible_site_str = ",".join(ai_possible_site_list)
                # 如果可能部位不匹配，直接跳过
                if not re.search(reg_site, ai_possible_site_str):
                    continue
                for ai_site in ai_site_list:
                    r = re.search("\D*(\d+)cm", ai_site, flags=re.I)
                    if r:
                        lenth = int(r.group(1))
                        if self.site_len_dict[reg_site][0] <= lenth <= self.site_len_dict[reg_site][1]:
                            self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                            matched_ai_index.append(j)
                            break

        # 第四次匹配 镜下所见中是部位，病理中是距离
        for reg_site in reg_site_list:
            if len(self.reg_site_ai_keys_dict[reg_site]) != 0:
                continue
            for j in range(len(self.classifed_ai_keys)):
                if j in matched_ai_index:
                    continue
                ai_obj_list = self.classifed_ai_keys[j]
                ai_site_list = [ai_obj.text for ai_obj in ai_obj_list if ai_obj.type == "site"]
                ai_possible_site_list = []
                for ai_obj in ai_obj_list:
                    if ai_obj.type == "site":
                        ai_possible_site_list += ai_obj.possible_site
                ai_possible_site_str = ",".join(ai_possible_site_list)
                # 直肠添加可能的部位 乙状结肠
                if re.search("直肠", ai_possible_site_str):
                    ai_possible_site_str += "," + "乙状结肠"
                if re.search(reg_site, ai_possible_site_str):
                    self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                    matched_ai_index.append(j)
                else:
                    r_ = re.search("(\d+)", reg_site)
                    if r_:
                        lenth = int(r_.group(1))
                        reg_site_ = None
                        for site_, lenth_list in self.site_len_dict.items():
                            if lenth_list[0] <= lenth <= lenth_list[1]:
                                reg_site_ = site_
                                break
                        if reg_site_ and re.search(reg_site_, ",".join(ai_site_list)):
                            self.reg_site_ai_keys_dict[reg_site] += ai_obj_list
                            matched_ai_index.append(j)
        # 第五次匹配（温州中心情况：病理中含有：直肠1，直肠2的）
        for reg_site in reg_site_list:
            if len(self.reg_site_ai_keys_dict[reg_site]) != 0:
                continue
            reg_site_ = reg_site[:-2]
            if reg_site_ not in reg_site_list:
                continue
            ai_obj_list = self.reg_site_ai_keys_dict[reg_site_]
            number_count = 0
            size_count = 0
            for ai_obj in ai_obj_list:
                if ai_obj.type == "size":
                    size_count += 1
                elif ai_obj.type == "number":
                    number_count += 1
            if number_count >= 2 or size_count >= 2:
                number_count = 0
                size_count = 0
                r = 0
                while r < len(ai_obj_list):
                    ai_obj = ai_obj_list[r]
                    if ai_obj.type == "size":
                        if number_count >= 1:
                            r += 1
                            break
                        size_count += 1
                    elif ai_obj.type == "number":
                        number_count += 1
                    r += 1
                if r < len(ai_obj_list):
                    self.reg_site_ai_keys_dict[reg_site_] = ai_obj_list[:r]
                    self.reg_site_ai_keys_dict[reg_site] = ai_obj_list[r:]

        # 特殊处理
        # 如果最后一个对象是number，并且number对象数量超过两个，则去除最后一个对象
        # 目的是处理： xx部位有3个xxcm的息肉，取1个送检 => 3个
        for reg_site, ai_obj_list in self.reg_site_ai_keys_dict.items():
            number_count = 0
            size_count = 0
            tmp_list = []
            for i in range(len(ai_obj_list)):
                ai_obj = ai_obj_list[i]
                if ai_obj.type == "number":
                    number_count += 1
                    tmp_list.append(['number', i])
                elif ai_obj.type == 'size':
                    size_count += 1
                    tmp_list.append(['size', i])
            if number_count >= 2 and size_count >= 1 and tmp_list[-1][0] == "number":
                ai_obj_list.pop(tmp_list[-1][1])
            self.reg_site_ai_keys_dict[reg_site] = ai_obj_list

    def clear_no_sonjian(self):
        """
        # 清除有送病检但没写的
        :return:
        """
        self.reg_site_ai_obj_list_dict = {}  # 结果
        for site, ai_obj_list in self.reg_site_ai_keys_dict.items():
            # 找到size对象和number对象的数量
            size_count = 0
            number_count = 0
            for ai_obj in ai_obj_list:
                if ai_obj.type == "size":
                    size_count += 1
                if ai_obj.type == "number":
                    number_count += 1
            # 找到有没有送检字样，用于排除没送检的
            songjian = False
            songjian_index_list = []
            for i in range(len(ai_obj_list)):
                ai_obj = ai_obj_list[i]
                # 如果有均送检，那么就不处理
                if ai_obj.type == "songjian" and not re.search("均", ai_obj.text):
                    songjian = True
                    songjian_index_list.append(i)

            obj_list_tmp = []
            # 存在送检字样，并且含有多对数量和大小的，就只找离送检最近的那个大小和数量
            if songjian and size_count >= 2 and number_count >= 2:
                for index in songjian_index_list:
                    site_count = 0
                    number_count = 0
                    while index >= 0:
                        if index > 0 and ai_obj_list[index - 1].type == "songjian":
                            break
                        if site_count != 0 and number_count != 0:
                            if ai_obj_list[index].type == "size" or ai_obj_list[index].type == "number":
                                pass
                            elif ai_obj_list[index].type == "site":
                                obj_list_tmp.append(ai_obj_list[index])
                            else:
                                break
                        else:
                            obj_list_tmp.append(ai_obj_list[index])
                            if ai_obj_list[index].type == "size":
                                site_count += 1
                            elif ai_obj_list[index].type == "number":
                                number_count += 1
                        index -= 1
            # 存在送检字样，但是不是多对，则数量对象只添加最大的
            elif songjian:
                max_number = 0
                for ai_obj in ai_obj_list:
                    if ai_obj.type != "number":
                        obj_list_tmp.append(ai_obj)
                    else:
                        if ai_obj.to_count() > max_number:
                            max_number = ai_obj.to_count()
                            obj_list_tmp.append(ai_obj)
            else:
                obj_list_tmp = ai_obj_list
            obj_list_tmp.sort(key=lambda x: x.start)
            if len(obj_list_tmp) > 0 and songjian:
                while len(obj_list_tmp) > 0 and obj_list_tmp[-1].type != "songjian":
                    obj_list_tmp.pop(-1)

            self.reg_site_ai_obj_list_dict[site] = obj_list_tmp

    def trans_cal(self):
        """
        # 转化格式并计算
        :return:
        """
        self.trans_cal_result = []  # 结果
        reg_site_ai_obj_list_dict = self.reg_site_ai_obj_list_dict
        tmp = []
        for reg_site, ai_obj_list in reg_site_ai_obj_list_dict.items():
            number_list = []
            size_list = []
            for i in range(len(ai_obj_list)):
                ai_obj = ai_obj_list[i]
                if ai_obj.type == "number":
                    number_list.append(ai_obj.to_count())
                elif ai_obj.type == "size":
                    size_list.append(ai_obj.get_size())
            sum_number = sum(number_list) if sum(number_list) != 0 else 0
            max_size = max(size_list) if len(size_list) != 0 else 0
            tmpa = {
                "site": reg_site,
                "number": sum_number,
                "max_size": max_size,
                "pathology": self.pathology_obj.get_site_desc(reg_site)
            }
            tmp.append(tmpa)
        # 如果数量小于病理部位数，则将为0的+1
        site_count = len(self.pathology_obj.site_reg_list)
        count = 0
        for item in tmp:
            count += item['number']
        if site_count > count:
            for item in tmp:
                if item['number'] == 0:
                    item['number'] = 1
        self.trans_cal_result = tmp

    def classify_ruzu(self):
        """
        是否入组
        :return:
        """
        if not self.pathology or self.pathology == "nan":
            self.ruzu = "无病理"
        elif self.classify_ca():
            self.ruzu = "癌"
        elif self.new_organism():
            self.ruzu = "新生物"
        elif self.classify_waike():
            self.ruzu = '肠外科手术、ESD术'
        elif self.classify_xirou_qiechu():
            self.ruzu = "非息肉相关疾病"
        elif self.younian_xirou():
            self.ruzu = "错构瘤性息肉"
        elif self.zhunbei_bujia():
            self.ruzu = '肠道准备不佳'
        elif self.size_error():
            self.ruzu = "大小书写有误"
        elif self.no_size():
            self.ruzu = '息肉没有写大小'
        elif not self.classify_match():
            self.ruzu = "病理不匹配"
        else:
            self.ruzu = '入组'

    def classify_ca(self):
        """
        判断是否为癌
        :return:
        """
        matches = re.finditer("癌|恶性", self.pathology)
        for match in matches:
            start = 0 if match.span()[0] - 10 < 0 else match.span()[0] - 10
            if not re.search("未见", self.pathology[start:match.span()[1]]):
                return 1
        return 0

    def new_organism(self):
        matchs = re.finditer("新生物|巨大隆起", self.microscopic_desc)
        for match in matchs:
            start = match.span()[0] - 15 if match.span()[0] >= 10 else 0
            if not re.search("未见", self.microscopic_desc[start:match.span()[1]]):
                return 1
        return 0

    def classify_waike(self):
        """
        肠外科手术、ESD术
        :return:
        """
        if re.search("肠CA.{0,8}术后|吻合口|根治术|吻合钉",
                     self.microscopic_diagnosis + self.microscopic_desc + self.pathology, flags=re.I):
            return 1
        if re.search("术后", self.microscopic_diagnosis) or re.search("(可见|呈)术后改变|ESD术后", self.microscopic_desc):
            if not re.search("APC术|APC治疗术后|电凝电切术后|息肉钳除术后|息肉切除术后", self.microscopic_diagnosis):
                return 1
        return 0

    def classify_xirou_qiechu(self):
        """
        判断是否为 非息肉相关疾病
        :return:
        """
        if re.search("脂肪瘤|梭形细胞瘤|神经内分泌|神经鞘|淋巴瘤", self.pathology):
            return 1
        if not re.search("息肉|隆起|腺瘤", self.microscopic_diagnosis + self.microscopic_desc + self.pathology, flags=0):
            return 1
        find = False
        for reg_site in self.pathology_obj.site_reg_list:
            if re.search("息肉|腺瘤|锯齿状|癌|无蒂|上皮内瘤变|瘤样增生|异型增生", self.pathology_obj.get_site_desc(reg_site)):
                find = True
        if not find:
            return 1
        return 0

    def younian_xirou(self):
        if re.search("P(-)?J.{0,8}综合征|家族性.{0,5}息肉病|(色素斑|色斑|色素|增生性|锯齿状)息肉综合征", self.microscopic_diagnosis,
                     flags=re.I):
            return True
        for reg_site in self.pathology_obj.site_reg_list:
            bingli_str = self.pathology_obj.get_site_desc(reg_site)
            if re.search("幼年性息肉", bingli_str):
                return True
        return False

    def zhunbei_bujia(self):
        if re.search("肠道.{0,6}准备.{0,6}(不佳|欠佳|不合格)|大量粪水残渣|影响观察视野|影响(视野|观察)", self.microscopic_desc):
            return True
        return False

    def no_size(self):
        for reg_site, ai_obj_list in self.reg_site_ai_keys_dict.items():
            if len(ai_obj_list) == 0:
                continue
            bingli_desc = self.pathology_obj.get_site_desc(reg_site)
            if re.search("腺瘤|锯齿状", bingli_desc):
                find_size = False
                for ai_obj in ai_obj_list:
                    if ai_obj.type == "size":
                        find_size = True
                if not find_size:
                    return True
        return False

    def size_error(self):
        """
        大小是否书写错误
        :return: True error
        """
        danwei = False
        exit_size = False
        for reg_site, ai_obj_list in self.reg_site_ai_obj_list_dict.items():
            for ai_obj in ai_obj_list:
                if ai_obj.type == "size":
                    matchs = re.finditer("(\d+\.?\d?)", str(ai_obj.text))
                    for match in matchs:
                        if str(match.group(0)).startswith("0") and "." not in str(match.group(0)):
                            return True
                    exit_size = True
                    if re.search("cm|mm|㎝", ai_obj.text, flags=re.I):
                        danwei = True
        if exit_size and not danwei:
            return True
        return False

    def classify_match(self):
        """
        判断病理是否匹配
        :return: True 匹配
        """
        for reg_site, ai_keys_list in self.reg_site_ai_obj_list_dict.items():
            if len(ai_keys_list) > 0:
                for ai_obj in ai_keys_list:
                    if ai_obj.type == "site":
                        if (not re.search("胃|十二指肠|食.{0,2}管|门齿|球.{0,2}部", ai_obj.text)) or re.search(
                                "乙状结肠|直肠|横结肠|升结肠|降结肠|结肠|回盲|肝曲|脾曲|盲肠", ai_obj.text):
                            return True
        return False

    def classify_level(self):
        """
        # 定级
        :return:
        """
        if self.level_standard == "Japan":
            self.level = JapanLevel(self).level
        elif self.level_standard == "Europe":
            self.level = EuropeLevel(self).level
        elif self.level_standard == "Britain":
            self.level = BritainLevel(self).level
        elif self.level_standard == "America":
            self.level = AmericaLevel(self).level
        elif self.level_standard == "China":
            self.level = ChinaLevel(self).level
