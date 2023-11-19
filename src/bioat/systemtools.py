# from bioat.api.internet_speed import main as speedtest


class SystemTools:
    def __init__(self):
        # self.test_internet_speed = speedtest()
        pass

#     def __convert_size(self, size, mode='b') -> float:
#         """将单位为Byte的数字转为其它单位
#
#         :param size: int, 单位为Byte的数字
#         :param mode: str, 单位，b、k、m、g、t或者B、K、M、G、T
#         :return:
#         """
#         assert mode in list('bkmgtBKMGT')
#         mode = mode.lower()
#
#         if mode == 'b':
#             return float(size)
#         elif mode == 'k':
#             return size / 1024.0
#         elif mode == 'm':
#             return size / 1024.0 / 1024.0
#         elif mode == 'g':
#             return size / 1024.0 / 1024.0 / 1024.0
#         elif mode == 't':
#             return size / 1024.0 / 1024.0 / 1024.0 / 1024.0
#         else:
#             return None
#
#     def __get_size(self, path, mode='b'):
#         """获取路径占用大小.
#         通过递归函数获取路径占用大小
#         :param path: str, 目标路径
#         :param mode: str, 单位，b、k、m、g、t或者B、K、M、G、T
#         :return: str, 占用大小
#         """
#         assert mode in list('bkmgtBKMGT')
#         mode = mode.lower()
#
#         # 目标是否存在
#         if os.path.exists(path):
#             # 目标路径是一个文件
#             if os.path.isfile(path):
#                 # 无序遍历，直接返回即可
#                 byte_size = os.path.getsize(path)
#                 return self.__convert_size(size=byte_size, mode=mode)
#
#             # 目标路径是一个文件夹
#             byte_size = 0.0
#             files = os.listdir(path)
#
#             for file in files:
#                 cur_path = os.path.join(path, file)
#                 try:
#                     byte_size += self.__get_size(cur_path)
#                 except TypeError:
#                     continue
#
#             return self.__convert_size(size=byte_size, mode=mode)
#         else:
#             # print(f'{path} not exists!')
#             pass
#
#
# if __name__ == '__main__':
#     # 大路径
#     aim_paths = [
#         # '/lustre1/chengqiyi_pkuhpc/zhaohn/2.disk_occupy/1T/fastq', # for test
#         # '/lustre1/chengqiyi_pkuhpc/zhaohn/2.disk_occupy/1T.1/1T/fastq', # for test\
#         #
#         '/lustre1/chengqiyi_pkuhpc/zhaohn',
#         '/lustre1/chengqiyi_pkuhpc/JiangZhe',
#         '/lustre1/chengqiyi_pkuhpc/LiuCong',
#         '/lustre1/chengqiyi_pkuhpc/liuwenqing',
#         '/lustre1/chengqiyi_pkuhpc/lubo',
#         '/lustre1/chengqiyi_pkuhpc/shx',
#         '/lustre2/chengqiyi_pkuhpc/bds',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/bds',
#         '/lustre2/chengqiyi_pkuhpc/dyy',
#         '/lustre2/chengqiyi_pkuhpc/LK',
#         '/lustre2/chengqiyi_pkuhpc/log',
#         '/lustre2/chengqiyi_pkuhpc/LSY',
#         '/lustre2/chengqiyi_pkuhpc/RenweiLi',
#         '/lustre2/chengqiyi_pkuhpc/wuhao',
#         '/lustre2/chengqiyi_pkuhpc/xiongxs_tmp',
#         '/lustre2/chengqiyi_pkuhpc/zhangxt',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/zhangxt',
#         '/lustre2/chengqiyi_pkuhpc/zml',
#         '/lustre3/chengqiyi_pkuhpc/Alumni/mshq',
#         '/lustre3/chengqiyi_pkuhpc/Alumni/Shuxiaoting',
#         '/lustre3/chengqiyi_pkuhpc/Alumni/SJH',
#         '/lustre3/chengqiyi_pkuhpc/Alumni/WK',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/dlt',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/Leizhixin',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/LiuJiangle',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/lixy',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/lzc',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/test2',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/ydy',
#         '/lustre3/chengqiyi_pkuhpc/folder_for_learning/ZhuangYuan',
#     ]
#
#     out_path = '/lustre1/chengqiyi_pkuhpc/zhaohn/__diskusage.all_user.csv'
#
#     dt_size = {key: [] for key in aim_paths}
#
#     for aim_path in aim_paths:
#         file_sizes = []
#         files = os.listdir(aim_path)
#
#         # print("aim_path: ", aim_path)
#         for file in files:
#             # continue
#             cur_path = os.path.join(aim_path, file)
#             # print("cur_path: ", cur_path)
#             if os.path.isdir(cur_path):
#                 try:
#                     file_size = (file, get_size(cur_path, 't'))
#                     file_sizes.append(file_size)
#                     print(f"INFO:    Processed {file_size[0]}  {file_size[1]:.3f}TB")
#                 except Exception as err:  # 捕捉所有异常
#                     print(f"WARNING: {err} @ {cur_path}, Skip ...")
#                 except FileNotFoundError:
#                     print(f"FileNotFoundError @ {cur_path}")
#                     continue
#                 except Exception as err:
#                     print(f"Exception @ {cur_path}: {err}")
#
#         # print("file_sizes: ", file_sizes)
#         file_sizes.sort(key=lambda t: t[1], reverse=True)
#         dt_size[aim_path] = file_sizes
#     date_time = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
#
#     with open(out_path, 'wt') as f:
#         f.write(f'{date_time},by,Huanan Zhao,\n,,,\n')
#         f.write("USER,PATH,TOTAL_SIZE,UNIT\n")
#         # print(dt_size)
#
#         for key in dt_size.keys():
#             user = os.path.basename(key)
#             size = sum([item[1] for item in dt_size[key]])
#             f.write(f'{user},{key},{size:>3f},TB\n')
#
#     print(f'one loop done! @ {date_time}')
