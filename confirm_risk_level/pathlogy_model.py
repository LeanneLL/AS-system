"""
@File ：fetch_gene_part.py
@Author ：wu_leanne
@Date ：2022/9/30
@Desc : This module is used to process pathological information
@Email : wu_leanne@163.com
"""
import re


class PathlogyRegDesc:
    """
    pathological model
    """

    def __init__(self, pathology):
        self.pathology = pathology
        self.site_reg_list = []
        self.congfu_site_reg_list = []
        self.site_reg_site_str_dict = {}
        self.site_reg_site_desc_dict = {}
        self.deal_pathology_data()

    def find_desc(self, start, pathology):
        """
        Look for a description of the pathological site
        :param start: starting index
        :param pathology: pathology
        :return:
        """
        end = start + 1
        while end < len(pathology):
            if re.search("（|\(|“|{", pathology[end]) and not re.search("亚蒂|整个|有蒂|无蒂", pathology[end:end + 3]):
                end += 1
                break
            end += 1
        return pathology[start:end]

    def deal_pathology_data(self):
        """
        Processing pathology data
        :return:
        """
        tmp_site_str_site_desc_dict = {}
        tmp_site_part_site_str_dict = {}

        # Common site
        common_site_list = ['乙状结肠', '直肠', '横结肠', '升结肠', '降结肠', '结肠', '回盲', '肝曲', '脾曲', "盲肠", "肛门内侧", '肝区', '脾区', "直乙",
                            "回肠末", "降乙"]
        pathology_list = self.pathology.split("\n")
        pathology_pattern_list = []  # Regular list of pathological sites (final list of results)
        pathology_list_part = []  # A list of strings containing sites
        r1 = re.search("手术部位(.*?)", self.pathology)
        if r1:
            r = re.search("手术部位(.*)", self.pathology)
            if r:
                pathology_list_part.append(r.group(1))
                desc = self.find_desc(r.span()[1], self.pathology)
                tmp_site_str_site_desc_dict[r.group(1)] = desc
        r2 = re.search("标本类型(.*?)", self.pathology)
        if r2:
            r = re.search("标本类型(.*)", self.pathology)
            if r:
                pathology_list_part.append(r.group(1))
                desc = self.find_desc(r.span()[1], self.pathology)
                tmp_site_str_site_desc_dict[r.group(1)] = desc
        for pathology in pathology_list:
            matchs = re.finditer(r"[（(“{](.+?)[）)”}]", pathology)
            for match in matchs:
                pathology_list_part.append(match.group(1))
                desc = self.find_desc(match.span()[1], pathology)
                tmp_site_str_site_desc_dict[match.group(1)] = desc
        # Splits a string containing sites
        pathology_list_part_ = []
        for pathology_part in pathology_list_part:
            if "范围约" in pathology_part:
                continue
            tmp_list = re.split("，|、|及", pathology_part)
            for i in tmp_list:
                pathology_list_part_.append(i)
                tmp_site_str_site_desc_dict[i] = tmp_site_str_site_desc_dict[pathology_part]
                tmp_site_part_site_str_dict[i] = pathology_part
        # The split string generates the final regular sites result
        for pathology_part in pathology_list_part_:
            re_ = re.search("肝曲|脾曲", pathology_part)
            re2_ = re.search("肝区|脾区", pathology_part)
            re3_ = re.search("近回盲部", pathology_part)
            if re_:
                # 如果有肝曲和脾曲，那么正则部位应该是 升结肠近肝曲 => 升结肠|肝曲
                pathology_part_ = re.sub("近", "", pathology_part)
                re_1 = re.search("肝曲|脾曲", pathology_part_)
                if re_1:
                    common_site = pathology_part_[:re_1.span()[0]] + "|" + pathology_part_[re_1.span()[0]:]
                    if common_site.startswith("|"):
                        common_site = common_site[1:]
                    common_site = common_site.split("曲")[0] + "曲"
                    if common_site not in pathology_pattern_list:
                        pathology_pattern_list.append(common_site)
                    else:
                        common_site = common_site + ".?"
                        self.congfu_site_reg_list.append(common_site)
                        pathology_pattern_list.append(common_site)
                    if not re.search("腺瘤", tmp_site_str_site_desc_dict.get(common_site, "")):
                        tmp_site_str_site_desc_dict[common_site] = tmp_site_str_site_desc_dict[pathology_part]
                    tmp_site_part_site_str_dict[common_site] = tmp_site_part_site_str_dict[pathology_part]
            elif re2_:
                # 如果有肝区和脾区,那么正则部位应该是 升结肠近肝区 => 升结肠|肝区
                pathology_part_ = re.sub("近", "", pathology_part)
                re_1 = re.search("肝区|脾区", pathology_part_)
                if re_1:
                    common_site = pathology_part_[:re_1.span()[0]] + "|" + pathology_part_[re_1.span()[0]:]
                    if common_site.startswith("|"):
                        common_site = common_site[1:]
                    common_site = common_site.split("区")[0] + "区"
                    if common_site not in pathology_pattern_list:
                        pathology_pattern_list.append(common_site)
                    else:
                        common_site = common_site + ".?"
                        self.congfu_site_reg_list.append(common_site)
                        pathology_pattern_list.append(common_site)
                    if not re.search("腺瘤", tmp_site_str_site_desc_dict.get(common_site, "")):
                        tmp_site_str_site_desc_dict[common_site] = tmp_site_str_site_desc_dict[pathology_part]
                    tmp_site_part_site_str_dict[common_site] = tmp_site_part_site_str_dict[pathology_part]
            elif re3_:
                # 如果有肝区和脾区,那么正则部位应该是 升结肠近肝区 => 升结肠|肝区
                pathology_part_ = re.sub("近", "", pathology_part)
                re_1 = re.search("回盲部", pathology_part_)
                if re_1:
                    common_site = pathology_part_[:re_1.span()[0]] + "|" + pathology_part_[re_1.span()[0]:]
                    if common_site.startswith("|"):
                        common_site = common_site[1:]
                    if common_site not in pathology_pattern_list:
                        pathology_pattern_list.append(common_site)
                    else:
                        common_site = common_site + ".?"
                        self.congfu_site_reg_list.append(common_site)
                        pathology_pattern_list.append(common_site)
                    if not re.search("腺瘤", tmp_site_str_site_desc_dict.get(common_site, "")):
                        tmp_site_str_site_desc_dict[common_site] = tmp_site_str_site_desc_dict[pathology_part]
                    tmp_site_part_site_str_dict[common_site] = tmp_site_part_site_str_dict[pathology_part]
            else:
                r_ = re.search(r"近", pathology_part)
                # 如果部位字符串中含有近xx，就 升结肠近横结肠 => 升结肠
                pathology_part_ = pathology_part[:r_.span(0)[0]] if r_ else pathology_part
                r1 = re.search("\D*(\d+)\D?cm", pathology_part_)
                if r1:
                    # 如果部位形式是距肛门xxcm
                    matchs = re.finditer(r"\D*(\d+)\D?cm", pathology_part)
                    for match in matchs:
                        site_ = match.group(1) + ".{0,6}(cm|㎝)"
                        if site_ not in pathology_pattern_list:
                            pathology_pattern_list.append(site_)
                        else:
                            site_1 = site_ + ".?"
                            self.congfu_site_reg_list.append(site_1)
                            pathology_pattern_list.append(site_1)
                            tmp_site_str_site_desc_dict[site_1] = tmp_site_str_site_desc_dict[site_]
                            tmp_site_part_site_str_dict[site_1] = tmp_site_part_site_str_dict[site_]

                        if not re.search("腺瘤", tmp_site_str_site_desc_dict.get(match.group(1) + ".{0,6}(cm|㎝)", "")):
                            tmp_site_str_site_desc_dict[site_] = tmp_site_str_site_desc_dict[
                                pathology_part]
                        tmp_site_part_site_str_dict[site_] = tmp_site_part_site_str_dict[
                            pathology_part]
                else:
                    # 不是xxcm形式
                    for common_site in common_site_list:
                        if common_site != "结肠" and re.search(common_site, pathology_part):
                            if common_site not in pathology_pattern_list:
                                pathology_pattern_list.append(common_site)
                            else:
                                common_site = common_site + ".?"
                                self.congfu_site_reg_list.append(common_site)
                                pathology_pattern_list.append(common_site)
                            if not re.search("腺瘤", tmp_site_str_site_desc_dict.get(common_site, "")):
                                tmp_site_str_site_desc_dict[common_site] = tmp_site_str_site_desc_dict[pathology_part]
                            tmp_site_part_site_str_dict[common_site] = tmp_site_part_site_str_dict[pathology_part]
                        elif common_site == "结肠" and (not re.search("乙状结肠|横结肠|升结肠|降结肠", pathology_part)) and (
                                re.search("结肠", pathology_part)):
                            if common_site not in pathology_pattern_list:
                                pathology_pattern_list.append(common_site)
                            else:
                                common_site = common_site + ".?"
                                self.congfu_site_reg_list.append(common_site)
                                pathology_pattern_list.append(common_site)
                            if not re.search("腺瘤", tmp_site_str_site_desc_dict.get(common_site, "")):
                                tmp_site_str_site_desc_dict[common_site] = tmp_site_str_site_desc_dict[pathology_part]
                            tmp_site_part_site_str_dict[common_site] = tmp_site_part_site_str_dict[pathology_part]
        # 如果按照上述方式没有找到部位，那么直接匹配
        if len(pathology_pattern_list) == 0:
            pat = "乙状结肠|直肠|横结肠|升结肠|降结肠|结肠|回盲|肝曲|脾曲|盲肠|肝区|脾区"
            matchs = re.finditer(pat, self.pathology)
            site_pat_list = []
            for match in matchs:
                site_pat_list.append({
                    "text": match.group(0),
                    "start": match.span()[0],
                    "end": match.span()[1]
                })
            for i in range(len(site_pat_list)):
                site_pat = site_pat_list[i]
                start = site_pat['end']
                if i < len(site_pat_list) - 1:
                    end = site_pat_list[i + 1]['start']
                else:
                    end = len(self.pathology)
                common_site = site_pat['text']
                if common_site not in pathology_pattern_list:
                    pathology_pattern_list.append(common_site)
                else:
                    common_site = common_site + ".?"
                    self.congfu_site_reg_list.append(common_site)
                    pathology_pattern_list.append(common_site)
                tmp_site_str_site_desc_dict[common_site] = self.pathology[start:end]
                tmp_site_part_site_str_dict[common_site] = site_pat['text']
        if len(pathology_pattern_list) == 0:
            matchs = re.finditer(r"\D*(\d+)\D?cm", self.pathology, flags=re.I)
            for match in matchs:
                common_site = match.group(1) + ".{0,6}(cm|㎝)"

                end = match.span()[1]
                while end <= len(self.pathology) - 1:
                    if self.pathology[end] == "\n":
                        break
                    end += 1
                start = match.span()[0]
                while start >= 0:
                    if self.pathology[start] == "\n":
                        break
                    start -= 1
                if common_site not in pathology_pattern_list:
                    pathology_pattern_list.append(common_site)
                else:
                    common_site = common_site + ".?"
                    self.congfu_site_reg_list.append(common_site)
                    pathology_pattern_list.append(common_site)
                tmp_site_str_site_desc_dict[common_site] = self.pathology[match.span()[0]:end]
                tmp_site_part_site_str_dict[common_site] = self.pathology[start:match.span()[0] + 1]
        # 区和曲都加入
        for i in range(len(pathology_pattern_list)):
            reg_site = pathology_pattern_list[i]
            if re.search("肝区", reg_site):
                pathology_pattern_list[i] = pathology_pattern_list[i] + "|肝曲"
            elif re.search("肝曲", reg_site):
                pathology_pattern_list[i] = pathology_pattern_list[i] + "|肝区"
            elif re.search("脾区", reg_site):
                pathology_pattern_list[i] = pathology_pattern_list[i] + "|脾曲"
            elif re.search("脾曲", reg_site):
                pathology_pattern_list[i] = pathology_pattern_list[i] + "|脾区"
            tmp_site_part_site_str_dict[pathology_pattern_list[i]] = tmp_site_part_site_str_dict[reg_site]
            tmp_site_str_site_desc_dict[pathology_pattern_list[i]] = tmp_site_str_site_desc_dict[reg_site]
        self.site_reg_list = list(set(pathology_pattern_list))  # 去重
        self.site_reg_list.sort(key=lambda x: [-len(x), x], reverse=True)  # 统一安装长度排序
        self.site_reg_site_str_dict = tmp_site_part_site_str_dict
        self.site_reg_site_desc_dict = tmp_site_str_site_desc_dict

    def get_site_str(self, site_reg):
        # 获取部位的原写法
        return self.site_reg_site_str_dict.get(site_reg)

    def get_site_desc(self, site_reg):
        # 获取部位描述
        return self.site_reg_site_desc_dict.get(site_reg)

    def __str__(self):
        return "|".join(self.site_reg_list)
