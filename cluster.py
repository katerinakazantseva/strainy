import csv
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from karateclub import LabelPropagation

bam="/Users/ekaterina.kazantseva/MT/test_data/test.bam"
snp="/Users/ekaterina.kazantseva/MT/test2.vcf"
edges=['edge_181','edge_11']
#edges=['edge_1','edge_2','edge_5', 'edge_598']
#edges=['edge_3', 'edge_4', 'edge_6', 'edge_7', 'edge_9', 'edge_10', 'edge_12', 'edge_13', 'edge_14', 'edge_15', 'edge_16', 'edge_17', 'edge_18', 'edge_19', 'edge_20', 'edge_21', 'edge_22', 'edge_23', 'edge_24', 'edge_25', 'edge_26', 'edge_27', 'edge_28', 'edge_29', 'edge_30', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_35', 'edge_36', 'edge_37', 'edge_38', 'edge_39', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_46', 'edge_47', 'edge_48', 'edge_49', 'edge_50', 'edge_51', 'edge_52', 'edge_53', 'edge_54', 'edge_55', 'edge_56', 'edge_57', 'edge_58', 'edge_59', 'edge_60', 'edge_61', 'edge_62', 'edge_63', 'edge_64', 'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_70', 'edge_71', 'edge_72', 'edge_73', 'edge_74', 'edge_75', 'edge_76', 'edge_77', 'edge_78', 'edge_79', 'edge_80', 'edge_81', 'edge_82', 'edge_83', 'edge_84', 'edge_85', 'edge_86', 'edge_87', 'edge_88', 'edge_89', 'edge_90', 'edge_91', 'edge_92', 'edge_93', 'edge_94', 'edge_95', 'edge_96', 'edge_97', 'edge_98', 'edge_99', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_104', 'edge_105', 'edge_106', 'edge_107', 'edge_108', 'edge_109', 'edge_110', 'edge_111', 'edge_112', 'edge_113', 'edge_114', 'edge_115', 'edge_116', 'edge_117', 'edge_118', 'edge_119', 'edge_120', 'edge_121', 'edge_122', 'edge_123', 'edge_124', 'edge_125', 'edge_126', 'edge_127', 'edge_128', 'edge_129', 'edge_130', 'edge_131', 'edge_132', 'edge_133', 'edge_134', 'edge_135', 'edge_136', 'edge_137', 'edge_138', 'edge_139', 'edge_140', 'edge_141', 'edge_142', 'edge_143', 'edge_144', 'edge_145', 'edge_146', 'edge_147', 'edge_148', 'edge_149', 'edge_150', 'edge_151', 'edge_152', 'edge_153', 'edge_154', 'edge_155', 'edge_156', 'edge_157', 'edge_158', 'edge_159', 'edge_160', 'edge_161', 'edge_162', 'edge_163', 'edge_164', 'edge_165', 'edge_166', 'edge_167', 'edge_168', 'edge_169', 'edge_170', 'edge_171', 'edge_172', 'edge_173', 'edge_174', 'edge_175', 'edge_176', 'edge_177', 'edge_178', 'edge_179', 'edge_180', 'edge_182', 'edge_183', 'edge_184', 'edge_185', 'edge_186', 'edge_187', 'edge_188', 'edge_189', 'edge_190', 'edge_191', 'edge_192', 'edge_193', 'edge_194', 'edge_195', 'edge_196', 'edge_197', 'edge_198', 'edge_199', 'edge_200', 'edge_201', 'edge_202', 'edge_203', 'edge_204', 'edge_205', 'edge_206', 'edge_207', 'edge_208', 'edge_209', 'edge_210', 'edge_211', 'edge_212', 'edge_213', 'edge_214', 'edge_215', 'edge_216', 'edge_217', 'edge_218', 'edge_219', 'edge_220', 'edge_221', 'edge_222', 'edge_223', 'edge_224', 'edge_225', 'edge_226', 'edge_227', 'edge_228', 'edge_229', 'edge_230', 'edge_231', 'edge_232', 'edge_233', 'edge_234', 'edge_235', 'edge_236', 'edge_237', 'edge_238', 'edge_239', 'edge_240', 'edge_241', 'edge_242', 'edge_243', 'edge_244', 'edge_245', 'edge_246', 'edge_247', 'edge_248', 'edge_249', 'edge_250', 'edge_251', 'edge_252', 'edge_253', 'edge_254', 'edge_255', 'edge_256', 'edge_257', 'edge_258', 'edge_259', 'edge_260', 'edge_261', 'edge_262', 'edge_263', 'edge_264', 'edge_265', 'edge_266', 'edge_267', 'edge_268', 'edge_269', 'edge_270', 'edge_271', 'edge_272', 'edge_273', 'edge_274', 'edge_275', 'edge_276', 'edge_277', 'edge_278', 'edge_279', 'edge_280', 'edge_281', 'edge_282', 'edge_283', 'edge_284', 'edge_285', 'edge_286', 'edge_287', 'edge_288', 'edge_289', 'edge_290', 'edge_291', 'edge_292', 'edge_293', 'edge_294', 'edge_295', 'edge_296', 'edge_297', 'edge_298', 'edge_299', 'edge_300', 'edge_301', 'edge_302', 'edge_303', 'edge_304', 'edge_305', 'edge_306', 'edge_307', 'edge_308', 'edge_309', 'edge_310', 'edge_311', 'edge_312', 'edge_313', 'edge_314', 'edge_315', 'edge_316', 'edge_317', 'edge_318', 'edge_319', 'edge_320', 'edge_321', 'edge_322', 'edge_323', 'edge_324', 'edge_325', 'edge_326', 'edge_327', 'edge_328', 'edge_329', 'edge_330', 'edge_331', 'edge_332', 'edge_333', 'edge_334', 'edge_335', 'edge_336', 'edge_337', 'edge_338', 'edge_339', 'edge_340', 'edge_341', 'edge_342', 'edge_343', 'edge_344', 'edge_345', 'edge_346', 'edge_347', 'edge_348', 'edge_349', 'edge_350', 'edge_351', 'edge_352', 'edge_353', 'edge_354', 'edge_355', 'edge_356', 'edge_357', 'edge_358', 'edge_359', 'edge_360', 'edge_361', 'edge_362', 'edge_363', 'edge_364', 'edge_365', 'edge_366', 'edge_367', 'edge_368', 'edge_369', 'edge_370', 'edge_371', 'edge_372', 'edge_373', 'edge_374', 'edge_375', 'edge_376', 'edge_377', 'edge_378', 'edge_379', 'edge_380', 'edge_381', 'edge_382', 'edge_383', 'edge_384', 'edge_385', 'edge_386', 'edge_387', 'edge_388', 'edge_389', 'edge_390', 'edge_391', 'edge_392', 'edge_393', 'edge_394', 'edge_395', 'edge_396', 'edge_397', 'edge_398', 'edge_399', 'edge_400', 'edge_401', 'edge_402', 'edge_403', 'edge_404', 'edge_405', 'edge_406', 'edge_407', 'edge_408', 'edge_409', 'edge_410', 'edge_411', 'edge_412', 'edge_413', 'edge_414', 'edge_415', 'edge_416', 'edge_417', 'edge_418', 'edge_419', 'edge_420', 'edge_421', 'edge_422', 'edge_423', 'edge_424', 'edge_425', 'edge_426', 'edge_427', 'edge_428', 'edge_429', 'edge_430', 'edge_431', 'edge_432', 'edge_433', 'edge_434', 'edge_435', 'edge_436', 'edge_437', 'edge_438', 'edge_439', 'edge_440', 'edge_441', 'edge_442', 'edge_443', 'edge_444', 'edge_445', 'edge_446', 'edge_447', 'edge_448', 'edge_449', 'edge_450', 'edge_451', 'edge_452', 'edge_453', 'edge_454', 'edge_455', 'edge_456', 'edge_457', 'edge_458', 'edge_459', 'edge_460', 'edge_461', 'edge_462', 'edge_463', 'edge_464', 'edge_465', 'edge_466', 'edge_467', 'edge_468', 'edge_469', 'edge_470', 'edge_471', 'edge_472', 'edge_473', 'edge_474', 'edge_475', 'edge_476', 'edge_477', 'edge_478', 'edge_479', 'edge_480', 'edge_481', 'edge_482', 'edge_483', 'edge_484', 'edge_485', 'edge_486', 'edge_487', 'edge_488', 'edge_489', 'edge_490', 'edge_491', 'edge_492', 'edge_493', 'edge_494', 'edge_495', 'edge_496', 'edge_497', 'edge_498', 'edge_499', 'edge_500', 'edge_501', 'edge_502', 'edge_503', 'edge_504', 'edge_505', 'edge_506', 'edge_507', 'edge_508', 'edge_509', 'edge_510', 'edge_511', 'edge_512', 'edge_513', 'edge_514', 'edge_515', 'edge_516', 'edge_517', 'edge_518', 'edge_519', 'edge_520', 'edge_521', 'edge_522', 'edge_523', 'edge_524', 'edge_525', 'edge_526', 'edge_527', 'edge_528', 'edge_529', 'edge_530', 'edge_531', 'edge_532', 'edge_533', 'edge_534', 'edge_535', 'edge_536', 'edge_537', 'edge_538', 'edge_539', 'edge_540', 'edge_541', 'edge_542', 'edge_543', 'edge_544', 'edge_545', 'edge_546', 'edge_547', 'edge_548', 'edge_549', 'edge_550', 'edge_551', 'edge_552', 'edge_553', 'edge_554', 'edge_555', 'edge_556', 'edge_557', 'edge_558', 'edge_559', 'edge_560', 'edge_561', 'edge_562', 'edge_563', 'edge_564', 'edge_565', 'edge_566', 'edge_567', 'edge_568', 'edge_569', 'edge_570', 'edge_571', 'edge_572', 'edge_573', 'edge_574', 'edge_575', 'edge_576', 'edge_577', 'edge_578', 'edge_579', 'edge_580', 'edge_581', 'edge_582', 'edge_583', 'edge_584', 'edge_585', 'edge_586', 'edge_587', 'edge_588', 'edge_589', 'edge_590', 'edge_591', 'edge_592', 'edge_593', 'edge_594', 'edge_595', 'edge_596', 'edge_597', 'edge_599', 'edge_600', 'edge_601', 'edge_602', 'edge_603', 'edge_604', 'edge_605', 'edge_606', 'edge_607', 'edge_608', 'edge_609', 'edge_610', 'edge_611', 'edge_612', 'edge_613', 'edge_614', 'edge_615', 'edge_616', 'edge_617', 'edge_618', 'edge_619', 'edge_620', 'edge_621', 'edge_622', 'edge_623', 'edge_624', 'edge_625', 'edge_626', 'edge_627', 'edge_628', 'edge_629', 'edge_630', 'edge_631', 'edge_632', 'edge_633', 'edge_634', 'edge_635', 'edge_636', 'edge_637', 'edge_638', 'edge_639', 'edge_640', 'edge_641', 'edge_642', 'edge_643', 'edge_644', 'edge_645', 'edge_646', 'edge_647', 'edge_648', 'edge_649', 'edge_650', 'edge_651', 'edge_652', 'edge_653', 'edge_654', 'edge_655', 'edge_656', 'edge_657', 'edge_658', 'edge_659', 'edge_660', 'edge_661', 'edge_662', 'edge_663', 'edge_664', 'edge_665', 'edge_666', 'edge_667', 'edge_668', 'edge_669', 'edge_670', 'edge_671', 'edge_672', 'edge_673', 'edge_674', 'edge_675', 'edge_676', 'edge_677', 'edge_678', 'edge_679', 'edge_680', 'edge_681', 'edge_682', 'edge_683', 'edge_684', 'edge_685', 'edge_686', 'edge_687', 'edge_688', 'edge_689', 'edge_690', 'edge_691', 'edge_692', 'edge_693', 'edge_694', 'edge_695', 'edge_696', 'edge_697', 'edge_698', 'edge_699', 'edge_700', 'edge_701', 'edge_702', 'edge_703', 'edge_704', 'edge_705', 'edge_706', 'edge_707', 'edge_708', 'edge_709', 'edge_710', 'edge_711', 'edge_712', 'edge_713', 'edge_714', 'edge_715', 'edge_716', 'edge_717', 'edge_718', 'edge_719', 'edge_720', 'edge_721', 'edge_722', 'edge_723', 'edge_724', 'edge_725', 'edge_726', 'edge_727', 'edge_728', 'edge_729', 'edge_730', 'edge_731', 'edge_732', 'edge_733', 'edge_734', 'edge_735', 'edge_736', 'edge_737', 'edge_738', 'edge_739', 'edge_740', 'edge_741', 'edge_742', 'edge_743', 'edge_744', 'edge_745', 'edge_746', 'edge_747', 'edge_748', 'edge_749', 'edge_750', 'edge_751', 'edge_752', 'edge_753', 'edge_754', 'edge_755', 'edge_756', 'edge_757', 'edge_758', 'edge_759', 'edge_760', 'edge_761', 'edge_762', 'edge_763', 'edge_764', 'edge_765', 'edge_766', 'edge_767', 'edge_768', 'edge_769', 'edge_770', 'edge_771', 'edge_772', 'edge_773', 'edge_774', 'edge_775', 'edge_776', 'edge_777', 'edge_778']
#edges=['edge_8']
#edges=['edge_57', 'edge_58', 'edge_59', 'edge_67', 'edge_68', 'edge_69', 'edge_70', 'edge_71', 'edge_72', 'edge_73', 'edge_74', 'edge_75', 'edge_85', 'edge_86', 'edge_89', 'edge_90', 'edge_91',  'edge_103',  'edge_105',  'edge_111',  'edge_118',  'edge_124',  'edge_131',  'edge_143',  'edge_145', 'edge_146', 'edge_147',  'edge_149', 'edge_150', 'edge_151',  'edge_161', 'edge_162', 'edge_163',  'edge_174', 'edge_175', 'edge_176',  'edge_178',  'edge_186', 'edge_187', 'edge_188', 'edge_189', 'edge_190', 'edge_191', 'edge_192', 'edge_193', 'edge_194', 'edge_195', 'edge_196', 'edge_197', 'edge_198', 'edge_199', 'edge_200', 'edge_201', 'edge_202', 'edge_203', 'edge_204', 'edge_205', 'edge_206', 'edge_207', 'edge_208', 'edge_209', 'edge_210', 'edge_211', 'edge_212', 'edge_213', 'edge_214', 'edge_215', 'edge_216', 'edge_217', 'edge_218', 'edge_219', 'edge_220', 'edge_221', 'edge_222', 'edge_223', 'edge_224', 'edge_225', 'edge_226', 'edge_227', 'edge_228', 'edge_229', 'edge_230', 'edge_231', 'edge_232', 'edge_233', 'edge_234', 'edge_235', 'edge_236', 'edge_237', 'edge_238', 'edge_239', 'edge_240', 'edge_241', 'edge_242', 'edge_243', 'edge_244', 'edge_245', 'edge_246', 'edge_247', 'edge_248', 'edge_249', 'edge_250', 'edge_251', 'edge_252', 'edge_253', 'edge_254', 'edge_255', 'edge_256', 'edge_257', 'edge_258', 'edge_259', 'edge_260', 'edge_261', 'edge_262', 'edge_263', 'edge_264', 'edge_265', 'edge_266', 'edge_267', 'edge_268', 'edge_269', 'edge_270', 'edge_271', 'edge_272', 'edge_273', 'edge_274', 'edge_275', 'edge_276', 'edge_277', 'edge_278', 'edge_279', 'edge_280', 'edge_281', 'edge_282', 'edge_283', 'edge_284', 'edge_285', 'edge_286', 'edge_287', 'edge_288', 'edge_289', 'edge_290', 'edge_291', 'edge_292', 'edge_293', 'edge_294', 'edge_295', 'edge_296', 'edge_297', 'edge_298', 'edge_299', 'edge_300', 'edge_301', 'edge_302', 'edge_303', 'edge_304', 'edge_305', 'edge_306', 'edge_307', 'edge_308', 'edge_309', 'edge_310', 'edge_311', 'edge_312', 'edge_313', 'edge_314', 'edge_315', 'edge_316', 'edge_317', 'edge_318', 'edge_319', 'edge_320', 'edge_321', 'edge_322', 'edge_323', 'edge_324', 'edge_325', 'edge_326', 'edge_327', 'edge_328', 'edge_329', 'edge_330', 'edge_331', 'edge_332', 'edge_333', 'edge_334', 'edge_335', 'edge_336', 'edge_337', 'edge_338', 'edge_339', 'edge_340', 'edge_341', 'edge_342', 'edge_343', 'edge_344', 'edge_345', 'edge_346', 'edge_347', 'edge_348', 'edge_349', 'edge_350', 'edge_351', 'edge_352', 'edge_353', 'edge_354', 'edge_355', 'edge_356', 'edge_357', 'edge_358', 'edge_359', 'edge_360', 'edge_361', 'edge_362', 'edge_363', 'edge_364', 'edge_365', 'edge_366', 'edge_367', 'edge_368', 'edge_369', 'edge_370', 'edge_371', 'edge_372', 'edge_373', 'edge_374', 'edge_375', 'edge_376', 'edge_377', 'edge_378', 'edge_379', 'edge_380', 'edge_381', 'edge_382', 'edge_383', 'edge_384', 'edge_385', 'edge_386', 'edge_387', 'edge_388', 'edge_389', 'edge_390', 'edge_391', 'edge_392', 'edge_393', 'edge_394', 'edge_395', 'edge_396', 'edge_397', 'edge_398', 'edge_399', 'edge_400', 'edge_401', 'edge_402', 'edge_403', 'edge_404', 'edge_405', 'edge_406', 'edge_407', 'edge_408', 'edge_409', 'edge_410', 'edge_411', 'edge_412', 'edge_413', 'edge_414', 'edge_415', 'edge_416', 'edge_417', 'edge_418', 'edge_419', 'edge_420', 'edge_421', 'edge_422', 'edge_423', 'edge_424', 'edge_425', 'edge_426', 'edge_427', 'edge_428', 'edge_429', 'edge_430', 'edge_431', 'edge_432', 'edge_433', 'edge_434', 'edge_435', 'edge_436', 'edge_437', 'edge_438', 'edge_439', 'edge_440', 'edge_441', 'edge_442', 'edge_443', 'edge_444', 'edge_445', 'edge_446', 'edge_447', 'edge_448', 'edge_449', 'edge_450', 'edge_451', 'edge_452', 'edge_453', 'edge_454', 'edge_455', 'edge_456', 'edge_457', 'edge_458', 'edge_459', 'edge_460', 'edge_461', 'edge_462', 'edge_463', 'edge_464', 'edge_465', 'edge_466', 'edge_467', 'edge_468', 'edge_469', 'edge_470', 'edge_471', 'edge_472', 'edge_473', 'edge_474', 'edge_475', 'edge_476', 'edge_477', 'edge_478', 'edge_479', 'edge_480', 'edge_481', 'edge_482', 'edge_483', 'edge_484', 'edge_485', 'edge_486', 'edge_487', 'edge_488', 'edge_489', 'edge_490', 'edge_491', 'edge_492', 'edge_493', 'edge_494', 'edge_495', 'edge_496', 'edge_497', 'edge_498', 'edge_499', 'edge_500', 'edge_501', 'edge_502', 'edge_503', 'edge_504', 'edge_505', 'edge_506', 'edge_507', 'edge_508', 'edge_509', 'edge_510', 'edge_511', 'edge_512', 'edge_513', 'edge_514', 'edge_515', 'edge_516', 'edge_517', 'edge_518', 'edge_519', 'edge_520', 'edge_521', 'edge_522', 'edge_523', 'edge_524', 'edge_525', 'edge_526', 'edge_527', 'edge_528', 'edge_529', 'edge_530', 'edge_531', 'edge_532', 'edge_533', 'edge_534', 'edge_535', 'edge_536', 'edge_537', 'edge_538', 'edge_539', 'edge_540', 'edge_541', 'edge_542', 'edge_543', 'edge_544', 'edge_545', 'edge_546', 'edge_547', 'edge_548', 'edge_549', 'edge_550', 'edge_551', 'edge_552', 'edge_553', 'edge_554', 'edge_555', 'edge_556', 'edge_557', 'edge_558', 'edge_559', 'edge_560', 'edge_561', 'edge_562', 'edge_563', 'edge_564', 'edge_565', 'edge_566', 'edge_567', 'edge_568', 'edge_569', 'edge_570', 'edge_571', 'edge_572', 'edge_573', 'edge_574', 'edge_575', 'edge_576', 'edge_577', 'edge_578', 'edge_579', 'edge_580', 'edge_581', 'edge_582', 'edge_583', 'edge_584', 'edge_585', 'edge_586', 'edge_587', 'edge_588', 'edge_589', 'edge_590', 'edge_591', 'edge_592', 'edge_593', 'edge_594', 'edge_595', 'edge_596', 'edge_597', 'edge_599', 'edge_600', 'edge_601', 'edge_602', 'edge_603', 'edge_604', 'edge_605', 'edge_606', 'edge_607', 'edge_608', 'edge_609', 'edge_610', 'edge_611', 'edge_612', 'edge_613', 'edge_614', 'edge_615', 'edge_616', 'edge_617', 'edge_618', 'edge_619', 'edge_620', 'edge_621', 'edge_622', 'edge_623', 'edge_624', 'edge_625', 'edge_626', 'edge_627', 'edge_628', 'edge_629', 'edge_630', 'edge_631', 'edge_632', 'edge_633', 'edge_634', 'edge_635', 'edge_636', 'edge_637', 'edge_638', 'edge_639', 'edge_640', 'edge_641', 'edge_642', 'edge_643', 'edge_644', 'edge_645', 'edge_646', 'edge_647', 'edge_648', 'edge_649', 'edge_650', 'edge_651', 'edge_652', 'edge_653', 'edge_654', 'edge_655', 'edge_656', 'edge_657', 'edge_658', 'edge_659', 'edge_660', 'edge_661', 'edge_662', 'edge_663', 'edge_664', 'edge_665', 'edge_666', 'edge_667', 'edge_668', 'edge_669', 'edge_670', 'edge_671', 'edge_672', 'edge_673', 'edge_674', 'edge_675', 'edge_676', 'edge_677', 'edge_678', 'edge_679', 'edge_680', 'edge_681', 'edge_682', 'edge_683', 'edge_684', 'edge_685', 'edge_686', 'edge_687', 'edge_688', 'edge_689', 'edge_690', 'edge_691', 'edge_692', 'edge_693', 'edge_694', 'edge_695', 'edge_696', 'edge_697', 'edge_698', 'edge_699', 'edge_700', 'edge_701', 'edge_702', 'edge_703', 'edge_704', 'edge_705', 'edge_706', 'edge_707', 'edge_708', 'edge_709', 'edge_710', 'edge_711', 'edge_712', 'edge_713', 'edge_714', 'edge_715', 'edge_716', 'edge_717', 'edge_718', 'edge_719', 'edge_720', 'edge_721', 'edge_722', 'edge_723', 'edge_724', 'edge_725', 'edge_726', 'edge_727', 'edge_728', 'edge_729', 'edge_730', 'edge_731', 'edge_732', 'edge_733', 'edge_734', 'edge_735', 'edge_736', 'edge_737', 'edge_738', 'edge_739', 'edge_740', 'edge_741', 'edge_742', 'edge_743', 'edge_744', 'edge_745', 'edge_746', 'edge_747', 'edge_748', 'edge_749', 'edge_750', 'edge_751', 'edge_752', 'edge_753', 'edge_754', 'edge_755', 'edge_756', 'edge_757', 'edge_758', 'edge_759', 'edge_760', 'edge_761', 'edge_762', 'edge_763', 'edge_764', 'edge_765', 'edge_766', 'edge_767', 'edge_768', 'edge_769', 'edge_770', 'edge_771', 'edge_772', 'edge_773', 'edge_774', 'edge_775', 'edge_776', 'edge_777', 'edge_778']

R=1





def read_snp(snp):
    SNP_pos = []
    vcf = open(snp, "rt")
    for line in vcf:
        if line.split()[0] == edge:
            SNP_pos.append(line.split()[1])

    #SNP_pos = ['492', '519', '533', '1287', '1373', '2746', '3346', '4027', '4531', '4597', '5125', '5149', '5164','5239', '5242', '5338', '5369', '7232', '7383', '8108', '8217']
    print(str(len(SNP_pos)) + " SNPs found")
    return(SNP_pos)


def read_bam(file,SNP_pos):
    bamfile = pysam.AlignmentFile(file, "rb")
    #iter = bamfile.fetch(edge)
    data = {}
    for pos in SNP_pos:
        for pileupcolumn in bamfile.pileup(edge, int(pos) - 1, int(pos), stepper='samtools', min_base_quality=0,
                                           ignore_overlaps=False,
                                           ignore_orphans=False, truncate=True):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                            pileupread.query_position]

                    except (KeyError):
                        data[pileupread.alignment.query_name] = {}
                        data[pileupread.alignment.query_name]["Start"] = pileupread.alignment.get_reference_positions()[
                            0]
                        data[pileupread.alignment.query_name]["Stop"] = pileupread.alignment.get_reference_positions()[
                            len(pileupread.alignment.get_reference_positions()) - 1]

                        data[pileupread.alignment.query_name][pos] = pileupread.alignment.query_sequence[
                            pileupread.query_position]

    bamfile.close()
    #print(data)
    return(data)

def distance(read1,read2,data, SNP_pos):
    d=-1
    for snp in SNP_pos:
        try:
            b1=data[read1][snp]
            b2=data[read2][snp]
            if b1 != b2 and len(b1)!=0 and  len(b2)!=0:
                if d==-1:
                    d=0
                d=d+1
            elif b1 == b2:
                if d==-1:
                    d=0
                d=d
        except:
            continue
        if d>=R:
            d=R
            break

        else:
            continue
    return (d)


def build_adj_matrix (cl,data,SNP_pos,I=1000):
    m = pd.DataFrame(-1, index=cl['ReadName'], columns=cl['ReadName'])
    for i in range(1,m.shape[1]):
        print(str(i)+"/"+str(m.shape[1])+" Reads processed \r" , end="")
        first_read=m.index[i]

        for j in range(0,i):
            second_read = m.index[j]

            #if len(set(data[first_read].keys()) & set(data[second_read].keys()))==0:
            if len(set(range(data[first_read]["Start"],data[first_read]["Stop"])) & set(range(data[second_read]["Start"],data[second_read]["Stop"])))<=I:
                m[second_read][first_read]=-1
            else:
                m[second_read][first_read] = distance(first_read,second_read, data,SNP_pos)
    return (m)


def remove_edges (m, R):
    m[m >= R] = -1
    return (m)

def change_w (m):
    m[m == 0] = 0.01
    m[m == -1] = 0
    m[m >=R] = 0
    return (m)


def main(edge):
    #READ SNPs
    print ("### Reading SNPs...")
    SNP_pos=read_snp(snp)

    #READ READS AND POSITIONS
    print ("### Reading Reads...")
    data=read_bam(bam,SNP_pos)
    cl=pd.DataFrame(data={'ReadName': data.keys()})
    cl['Cluster'] = 'NA'
    print (str(len(cl['ReadName']))+" reads found")


    #CALCULATE DISTANCE and ADJ MATRIX
    print ("### Calculatind distances/Building adj matrix...")
    m=build_adj_matrix (cl,data,SNP_pos)
    m.to_csv("output/adj_M_%s.csv" % edge)
    m=pd.read_csv("output/adj_M_%s.csv" % edge,index_col='ReadName')
    #print("### Removing overweighed egdes...")
    m=remove_edges (m, R)


    #BUILD graph
    print("### Creating graph...")
    m1=m
    #print(m)
    #print(cl)
    m1.columns = range(0,len(cl['ReadName']))
    m1.index=range(0,len(cl['ReadName']))
    G = nx.from_pandas_adjacency(change_w(m1.transpose()))
    to_remove = [(a, b) for a, b, attrs in G.edges(data=True) if  attrs["weight"] == 0]
    G.remove_edges_from(to_remove)
    nx.draw(G, with_labels = False, width=0.03,node_color='pink', node_size=3, font_size=5)
    ln=pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[4]
    cov=pysam.samtools.coverage("-r", edge, bam, "--no-header").split()[8]


    print("### Searching clusters...")


    #cliques=sorted(list(nx.find_cliques(G)),key=len, reverse=True)
    model = LabelPropagation()
    model.fit(G)
    cluster_membership = model.get_memberships()
    #print("CLUSTERS")
    clN = 0
    uncl = 0

    for value in set(cluster_membership.values()):
        group=[k for k,v in cluster_membership.items() if v == value]
        if len(group)>2:
            clN=clN+1
            #print("cluster " + str(value))
            #print(len(group))
            #print(cl['ReadName'][group])
            cl['Cluster'][group]=value
        else:
            uncl=uncl+1
    cl.to_csv("output/clusters_%s.csv" % edge)

    plt.suptitle(str(edge)+" coverage:"+str(cov)+" length:"+str(ln)+" clN:"+str(clN))
    #plt.savefig("output/graph_%s_I1000.png" % edge, format="PNG")
    plt.savefig("output/graph_%s_I1000.png" % edge, format="PNG")
    plt.close()
    #add tags to bam


    def write_bam(file, cl,edge):
        bamfile = pysam.AlignmentFile(file, "rb")
        for pileupcolumn in bamfile.pileup(edge):
            for pileupread in pileupcolumn.pileups:
                #clN=cl.loc[cl['ReadName'] == pileupread.alignment.query_name].values
                clN=1
                pileupread.alignment.tags+=[('cluster', clN)]
        infile = pysam.AlignmentFile("-", "rb")
        fo = pysam.Samfile('/output/test2edge8.bam', "wb",template=infile)
        fo.write(Read)


    #write_bam(bam,cl,edge)


    # Calculate statistics

    print("Summary for: " + edge)
    print("Clusters found: " + str(clN))
    print("Reads unclassified: " + str(uncl))
    print("Number of reads in each cluster: ")
    print(cl['Cluster'].value_counts(dropna=False))


for edge in edges:
    print("------------- "+str(edge)+" ------------- ")
    main(edge)















#UNUSED PART


print("### Searching connected components...")


def find(read_name, clN):

        x = m.loc[m[read_name].isin((0,R-1))].index
        if cl[(cl['ReadName'] == read_name)]['Cluster'].values==['NA']:
            #print("assing cl "+str(clN))
            #print(read_name)
            #print(x)
            cl.loc[cl['ReadName'] == read_name,'Cluster']=clN
            for i in x:
                find(i,clN)





#for next_read in cl['ReadName']:
 #   if cl[(cl['ReadName'] == next_read)]['Cluster'].values==['NA']:
  #      x = m.loc[m[next_read].isin((0,R-1))].index
   #     if len(x)>0:
    #        clN = clN + 1
     #       find(next_read,clN)

      #  else:
       #     uncl = uncl + 1
        #    print(str(next_read)+" unclustered")
         #   cl.loc[cl['ReadName'] == next_read, 'Cluster'] = 'unclustered'




#cl.to_csv("clusters-v2_%s.csv" % edge)
#cl=pd.read_csv("clusters-v2_edge_8.csv")
#print(df)

















