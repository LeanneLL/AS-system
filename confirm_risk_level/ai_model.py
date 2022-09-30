"""
@File ：fetch_gene_part.py
@Author ：wu_leanne
@Date ：2022/9/30
@Desc : This module is used to encapsulate AS model results as objects
@Email : wu_leanne@163.com
"""

import re
import cn2an


class AiObject:
    """
    base model
    """

    def __init__(self, text="", start=0, end=0):
        """
        :param text: microscopic description message
        :param start: starting index
        :param end: ending index
        """
        self.text = text
        self.start = start
        self.end = end

    def to_count(self):
        """
        The number of polyps was obtained by converting Chinese numbers into Arabic numbers
        :return: Number of polyps
        """
        if self.text == "密集":
            return 10
        if re.search("数枚|多枚|数个|多发|多个|数颗|多处|散在", self.text):
            return 5
        cn_num = str(self.text)
        if cn_num.isdigit():
            return int(cn_num)
        cn_num = cn2an.transform(cn_num)
        r_ = re.search("(\d+)", cn_num)
        if r_:
            return int(r_.group(0))
        return 0

    def get_size(self):
        """
        The maximum polyp size was obtained
        :return: maximum size
        """
        size = 0
        mm = re.search("mm", str(self.text), flags=re.I)
        matchs = re.finditer("(\d+\.?\d?)", str(self.text))
        for match in matchs:
            tmp_size = float(match.group(0)) / 10 if mm else float(match.group(0))
            if str(match.group(0)).startswith("0") and "." not in str(match.group(0)):
                tmp_size = tmp_size / 10
            size = max(size, tmp_size)
        return size


class AiSite(AiObject):
    """
    site model
    """

    def __init__(self, text="", start=0, end=0):
        super().__init__(text, start, end)
        self.type = "site"
        self.possible_site = []

    def __str__(self):
        return "{" + f"'text': '{self.text}', 'start': {self.start}, 'end': {self.end}, 'type': '{self.type}', 'possible_site': '{self.possible_site}'" + "}"


class AiNumber(AiObject):
    """
    count model
    """

    def __init__(self, text="", start=0, end=0):
        super().__init__(text, start, end)
        self.type = "number"

    def __str__(self):
        return "{" + f"'text': '{self.text}', 'start': {self.start}, 'end': {self.end}, 'type': '{self.type}'" + "}"


class AiSize(AiObject):
    """
    size model
    """

    def __init__(self, text="", start=0, end=0):
        super().__init__(text, start, end)
        self.type = "size"

    def __str__(self):
        return "{" + f"'text': '{self.text}', 'start': {self.start}, 'end': {self.end}, 'type': '{self.type}'" + "}"


class AiBiopsy(AiObject):
    """
    biopsy model
    """

    def __init__(self, text="", start=0, end=0):
        super().__init__(text, start, end)
        self.type = "songjian"

    def __str__(self):
        return "{" + f"'text': '{self.text}', 'start': {self.start}, 'end': {self.end}, 'type': '{self.type}'" + "}"
