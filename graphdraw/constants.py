DPI = 400
REF_LEN = {'chr1': 248387328,
           'chr2': 242696752,
           'chr3': 201105948,
           'chr4': 193574945,
           'chr5': 182045439,
           'chr6': 172126628,
           'chr7': 160567428,
           'chr8': 146259331,
           'chr9': 150617247,
           'chr10': 134758134,
           'chr11': 135127769,
           'chr12': 133324548,
           'chr13': 113566686,
           'chr14': 101161492,
           'chr15': 99753195,
           'chr16': 96330374,
           'chr17': 84276897,
           'chr18': 80542538,
           'chr19': 61707364,
           'chr20': 66210255,
           'chr21': 45090682,
           'chr22': 51324926,
           'chrX': 154259566,
           'chrY': 62460029,
           'chrM': 16569}


"""
The centromere locations for t2t, hg19, and hg38 were taken from the Chromosome band (ideogram) tables from ucsc table browser, and then the locations with acen tag were taken.

hg38 has another table called centromeres but has several entries for each chromosome and was not able to figure out the differences, I also found other tables with centromeric locations, but I will stick to the cytoband table for consistency.

As for gaps, they were also taken table browser using the GAP table and took all entries
"""

CENTRO_T2T_UCSC_IDEOGRAM = {'chr1': (121796048, 126300487),
                            'chr2': (92333543, 94673023),
                            'chr3': (91738002, 96415026),
                            'chr4': (49705154, 55199795),
                            'chr5': (47039134, 49596625),
                            'chr6': (58286706, 61058390),
                            'chr7': (60414372, 63714499),
                            'chr8': (44215832, 46325080),
                            'chr9': (44951775, 47582595),
                            'chr10': (39633793, 41664589),
                            'chr11': (51035789, 54450838),
                            'chr12': (34620838, 37202490),
                            'chr13': (15547593, 17498291),
                            'chr14': (10092112, 12708411),
                            'chr15': (16678794, 17694466),
                            'chr16': (35848286, 37829521),
                            'chr17': (23892419, 27486939),
                            'chr18': (15965699, 20933550),
                            'chr19': (25817676, 29768171),
                            'chr20': (26925852, 29099655),
                            'chr21': (10962853, 11306205),
                            'chr22': (12788180, 15711065),
                            'chrX': (57820107, 60927025),
                            'chrM': (0, 0),
                            'chrY': (10565750, 10883085)}

CENTRO_HG19_UCSC_IDEOGRAM = {'chr1': (121500000, 128900000),
                             'chr2': (90500000, 96800000),
                             'chr3': (87900000, 93900000),
                             'chr4': (48200000, 52700000),
                             'chr5': (46100000, 50700000),
                             'chr6': (58700000, 63300000),
                             'chr7': (58000000, 61700000),
                             'chr8': (43100000, 48100000),
                             'chr9': (47300000, 50700000),
                             'chrX': (58100000, 63000000),
                             'chrY': (11600000, 13400000),
                             'chr10': (38000000, 42300000),
                             'chr11': (51600000, 55700000),
                             'chr12': (33300000, 38200000),
                             'chr13': (16300000, 19500000),
                             'chr14': (16100000, 19100000),
                             'chr15': (15800000, 20700000),
                             'chr16': (34600000, 38600000),
                             'chr17': (22200000, 25800000),
                             'chr18': (15400000, 19000000),
                             'chr19': (24400000, 28600000),
                             'chr21': (10900000, 14300000),
                             'chr22': (12200000, 17900000)}


CENTRO_HG38_UCSC_IDEOGRAM = {'chr1': (121700000, 125100000),
                             'chr2': (91800000, 96000000),
                             'chr3': (87800000, 94000000),
                             'chr4': (48200000, 51800000),
                             'chr5': (46100000, 51400000),
                             'chr6': (58500000, 62600000),
                             'chr7': (58100000, 62100000),
                             'chr8': (43200000, 47200000),
                             'chr9': (42200000, 45500000),
                             'chrX': (58100000, 63800000),
                             'chrY': (10300000, 10600000),
                             'chr10': (38000000, 41600000),
                             'chr11': (51000000, 55800000),
                             'chr12': (33200000, 37800000),
                             'chr13': (16500000, 18900000),
                             'chr14': (16100000, 18200000),
                             'chr15': (17500000, 20500000),
                             'chr16': (35300000, 38400000),
                             'chr17': (22700000, 27400000),
                             'chr18': (15400000, 21500000),
                             'chr19': (24200000, 28100000),
                             'chr20': (25700000, 30400000),
                             'chr21': (10900000, 13000000),
                             'chr22': (13700000, 17400000)}


HG19_GAPS = {'chr1': [(0, 10000),
  (177417, 227417),
  (267719, 317719),
  (471368, 521368),
  (2634220, 2684220),
  (3845268, 3995268),
  (13052998, 13102998),
  (13219912, 13319912),
  (13557162, 13607162),
  (17125658, 17175658),
  (29878082, 30028082),
  (103863906, 103913906),
  (120697156, 120747156),
  (120936695, 121086695),
  (121485434, 121535434),
  (121535434, 124535434),
  (124535434, 142535434),
  (142731022, 142781022),
  (142967761, 143117761),
  (143292816, 143342816),
  (143544525, 143644525),
  (143771002, 143871002),
  (144095783, 144145783),
  (144224481, 144274481),
  (144401744, 144451744),
  (144622413, 144672413),
  (144710724, 144810724),
  (145833118, 145883118),
  (146164650, 146214650),
  (146253299, 146303299),
  (148026038, 148176038),
  (148361358, 148511358),
  (148684147, 148734147),
  (148954460, 149004460),
  (149459645, 149509645),
  (205922707, 206072707),
  (206332221, 206482221),
  (223747846, 223797846),
  (235192211, 235242211),
  (248908210, 249058210),
  (249240621, 249250621)],
 'chr10': [(0, 10000),
  (10000, 60000),
  (17974675, 18024675),
  (38818835, 38868835),
  (39154935, 39254935),
  (39254935, 42254935),
  (42254935, 42354935),
  (42546687, 42596687),
  (46426964, 46476964),
  (47429169, 47529169),
  (47792476, 47892476),
  (48055707, 48105707),
  (49095536, 49195536),
  (51137410, 51187410),
  (51398845, 51448845),
  (125869472, 125919472),
  (128616069, 128766069),
  (133381404, 133431404),
  (133677527, 133727527),
  (135524747, 135534747)],
 'chr10_jh591181_fix': [(2010264, 2060264)],
 'chr11': [(0, 10000),
  (10000, 60000),
  (1162759, 1212759),
  (50783853, 50833853),
  (50833853, 51040853),
  (51040853, 51090853),
  (51594205, 51644205),
  (51644205, 54644205),
  (54644205, 54694205),
  (69089801, 69139801),
  (69724695, 69774695),
  (87688378, 87738378),
  (96287584, 96437584),
  (134946516, 134996516),
  (134996516, 135006516)],
 'chr11_jh806581_fix': [(499909, 549909)],
 'chr12': [(0, 10000),
  (10000, 60000),
  (95739, 145739),
  (7189876, 7239876),
  (34856694, 37856694),
  (109373470, 109423470),
  (122530623, 122580623),
  (132706992, 132806992),
  (133841895, 133851895)],
 'chr13': [(0, 10000),
  (10000, 16000000),
  (16000000, 19000000),
  (19000000, 19020000),
  (86760324, 86910324),
  (112353994, 112503994),
  (114325993, 114425993),
  (114639948, 114739948),
  (115109878, 115159878),
  (115159878, 115169878)],
 'chr14': [(0, 10000),
  (10000, 16000000),
  (16000000, 19000000),
  (107289540, 107339540),
  (107339540, 107349540)],
 'chr15': [(0, 10000),
  (10000, 17000000),
  (17000000, 20000000),
  (20894633, 20935075),
  (21398819, 21885000),
  (22212114, 22262114),
  (22596193, 22646193),
  (23514853, 23564853),
  (29159443, 29209443),
  (82829645, 82879645),
  (84984473, 85034473),
  (102521392, 102531392)],
 'chr16': [(0, 10000),
  (10000, 60000),
  (8636921, 8686921),
  (34023150, 34173150),
  (35285801, 35335801),
  (35335801, 38335801),
  (38335801, 46335801),
  (46335801, 46385801),
  (88389383, 88439383),
  (90294753, 90344753),
  (90344753, 90354753)],
 'chr17': [(296626, 396626),
  (21566608, 21666608),
  (22263006, 25263006),
  (34675848, 34725848),
  (62410760, 62460760),
  (77546461, 77596461),
  (79709049, 79759049)],
 'chr17_ctg5_hap1': [(1256794, 1306794), (1588968, 1638968)],
 'chr17_gl383561_fix': [(448214, 466467),
  (513301, 513321),
  (628875, 628895),
  (635816, 635918),
  (636913, 637464)],
 'chr17_jh720447_fix': [(190775, 240775)],
 'chr17_ke332502_fix': [(100012, 100085),
  (101792, 102078),
  (104926, 105301),
  (160048, 160456),
  (181878, 182834),
  (196662, 197232),
  (267213, 267233),
  (281523, 281724),
  (321182, 321624),
  (338198, 338308)],
 'chr18': [(0, 10000),
  (15410898, 15460898),
  (15460898, 18460898),
  (18460898, 18510898),
  (52059136, 52209136),
  (72283353, 72333353),
  (75721820, 75771820),
  (78017248, 78067248),
  (78067248, 78077248)],
 'chr18_gl383571_alt': [(44225, 94225)],
 'chr19': [(0, 10000),
  (10000, 60000),
  (7346004, 7396004),
  (8687198, 8737198),
  (20523415, 20573415),
  (24631782, 24681782),
  (24681782, 27681782),
  (27681782, 27731782),
  (59118983, 59128983)],
 'chr19_gl949746_alt': [(309793, 314793), (700652, 800652)],
 'chr19_gl949747_alt': [(195631, 256631), (394094, 399094)],
 'chr19_gl949748_alt': [(195631, 720631), (761613, 811613), (845239, 877239)],
 'chr19_gl949749_alt': [(195631, 807631), (860776, 904776)],
 'chr19_gl949750_alt': [(195631, 819631), (855325, 879325)],
 'chr19_gl949751_alt': [(195631, 705631), (789618, 815618)],
 'chr19_gl949752_alt': [(253362, 258362)],
 'chr19_gl949753_alt': [(195631, 213631), (361077, 363077), (509414, 609414)],
 'chr2': [(0, 10000),
  (3529312, 3579312),
  (5018788, 5118788),
  (16279724, 16329724),
  (21153113, 21178113),
  (87668206, 87718206),
  (89630436, 89830436),
  (90321525, 90371525),
  (90545103, 91545103),
  (91545103, 91595103),
  (92326171, 95326171),
  (110109337, 110251337),
  (149690582, 149790582),
  (234003741, 234053741),
  (239801978, 239831978),
  (240784132, 240809132),
  (243102476, 243152476),
  (243189373, 243199373)],
 'chr20': [(0, 10000),
  (10000, 60000),
  (26319569, 26369569),
  (26369569, 29369569),
  (29369569, 29419569),
  (29653908, 29803908),
  (34897085, 34947085),
  (61091437, 61141437),
  (61213369, 61263369),
  (62965520, 63015520),
  (63015520, 63025520)],
 'chr21': [(0, 10000),
  (10000, 5211193),
  (5211193, 9411193),
  (9595548, 9645548),
  (9775437, 9825437),
  (10034920, 10084920),
  (10215976, 10365976),
  (10647896, 10697896),
  (11188129, 11238129),
  (11238129, 11288129),
  (11288129, 14288129),
  (14288129, 14338129),
  (42955559, 43005559),
  (44632664, 44682664),
  (48119895, 48129895)],
 'chr22': [(0, 10000),
  (10000, 13000000),
  (13000000, 16000000),
  (16000000, 16050000),
  (16697850, 16847850),
  (20509431, 20609431),
  (50364777, 50414777),
  (51244566, 51294566),
  (51294566, 51304566)],
 'chr3': [(0, 10000),
  (10000, 60000),
  (66170270, 66270270),
  (90504854, 93504854),
  (194041961, 194047251),
  (197962430, 198012430),
  (198012430, 198022430)],
 'chr4': [(0, 10000),
  (1423146, 1478646),
  (8799203, 8818203),
  (9274642, 9324642),
  (31820917, 31837417),
  (32834638, 32840638),
  (40296396, 40297096),
  (49338941, 49488941),
  (49660117, 52660117),
  (59739333, 59789333),
  (75427379, 75452279),
  (191044276, 191144276),
  (191144276, 191154276)],
 'chr4_gl877872_fix': [(216051, 266051)],
 'chr5': [(0, 10000),
  (17530657, 17580657),
  (46405641, 49405641),
  (91636128, 91686128),
  (138787073, 138837073),
  (155138727, 155188727),
  (180905260, 180915260)],
 'chr5_gl339449_alt': [(1138692, 1238692)],
 'chr6': [(0, 10000),
  (10000, 60000),
  (58087659, 58137659),
  (58780166, 58830166),
  (58830166, 61830166),
  (61830166, 61880166),
  (62128589, 62178589),
  (95680543, 95830543),
  (157559467, 157609467),
  (157641300, 157691300),
  (167942073, 168042073),
  (170279972, 170329972),
  (171055067, 171105067),
  (171105067, 171115067)],
 'chr6_apd_hap1': [(207807, 252060),
  (327074, 356458),
  (448574, 536172),
  (588761, 692953),
  (944435, 963558),
  (970040, 1029417),
  (1173934, 1198528),
  (1307951, 1317697),
  (1341208, 1344139),
  (1451581, 1556172),
  (1582549, 1598989),
  (1736596, 1786349),
  (1833811, 1937682),
  (2009922, 2164480),
  (2169002, 2254169),
  (2359974, 2761234),
  (2821679, 2858431),
  (2894657, 3037482),
  (3122973, 3136480),
  (3180323, 3355295),
  (3370179, 3415446),
  (3491066, 3594101),
  (3638667, 3642320),
  (3696623, 3736194),
  (3762870, 3766192),
  (3790366, 3970602),
  (4220467, 4274176),
  (4390029, 4597885)],
 'chr6_dbb_hap3': [(138055, 245960),
  (622457, 640380),
  (1337861, 1344119),
  (1836690, 1837200),
  (2117295, 2119813),
  (2606257, 2723615),
  (3407000, 3481304),
  (3736386, 3772614),
  (4206542, 4226964),
  (4354126, 4376794)],
 'chr6_mann_hap4': [(145852, 184243),
  (1388812, 1438812),
  (2047420, 2056121),
  (2252213, 2277026),
  (2312291, 2324610),
  (2517052, 2540962),
  (2772006, 2810953),
  (2921096, 2929556),
  (3062107, 3151453),
  (3174271, 3227203),
  (3261361, 3374492),
  (3412258, 3450867),
  (3553455, 3557741),
  (3895092, 3926341),
  (4236950, 4245095),
  (4441658, 4480941)],
 'chr6_mcf_hap5': [(135141, 138043),
  (179429, 186757),
  (229133, 282732),
  (383988, 404415),
  (491590, 560421),
  (676599, 685588),
  (994250, 1009296),
  (1177759, 1292220),
  (1562952, 1726699),
  (1815769, 1838106),
  (2091199, 2189158),
  (2241824, 2251998),
  (2300850, 2305647),
  (2712451, 2811903),
  (2947576, 2958999),
  (3154498, 3188719),
  (3329672, 3355524),
  (3602941, 3659279),
  (3899294, 3929849),
  (4149579, 4149626),
  (4284353, 4302274),
  (4422172, 4594253)],
 'chr6_qbl_hap6': [(521900, 681210),
  (2680914, 2732064),
  (3023176, 3049607),
  (3316246, 3369039),
  (3993031, 4020006)],
 'chr6_ssto_hap7': [(861800, 941242),
  (1727717, 1744408),
  (1836632, 1838892),
  (1983561, 1997824),
  (2052735, 2087415),
  (2101474, 2113875),
  (2138491, 2303018),
  (2411247, 2461193),
  (3040184, 3054966),
  (3093568, 3154699),
  (3191196, 3214366),
  (3368931, 3445987),
  (4370564, 4420564),
  (4580410, 4610218),
  (4639806, 4661556),
  (4723509, 4774367),
  (4783798, 4836049)],
 'chr7': [(0, 10000),
  (232484, 282484),
  (50370631, 50410631),
  (58054331, 61054331),
  (61310513, 61360513),
  (61460465, 61510465),
  (61677020, 61727020),
  (61917157, 61967157),
  (74715724, 74765724),
  (100556043, 100606043),
  (130154523, 130254523),
  (139379377, 139404377),
  (142048195, 142098195),
  (142276197, 142326197),
  (143347897, 143397897),
  (154270634, 154370634),
  (159128663, 159138663)],
 'chr7_jh636058_fix': [(170100, 220100)],
 'chr8': [(0, 10000),
  (7474649, 7524649),
  (12091854, 12141854),
  (43838887, 46838887),
  (48130499, 48135599),
  (86576451, 86726451),
  (142766515, 142816515),
  (145332588, 145432588),
  (146304022, 146354022),
  (146354022, 146364022)],
 'chr9': [(0, 10000),
  (39663686, 39713686),
  (39974796, 40024796),
  (40233029, 40283029),
  (40425834, 40475834),
  (40940341, 40990341),
  (41143214, 41193214),
  (41365793, 41415793),
  (42613955, 42663955),
  (43213698, 43313698),
  (43946569, 43996569),
  (44676646, 44726646),
  (44908293, 44958293),
  (45250203, 45350203),
  (45815521, 45865521),
  (46216430, 46266430),
  (46461039, 46561039),
  (47060133, 47160133),
  (47317679, 47367679),
  (47367679, 50367679),
  (50367679, 65367679),
  (65367679, 65467679),
  (65918360, 65968360),
  (66192215, 66242215),
  (66404656, 66454656),
  (66614195, 66664195),
  (66863343, 66913343),
  (67107834, 67207834),
  (67366296, 67516296),
  (67987998, 68137998),
  (68514181, 68664181),
  (68838946, 68988946),
  (69278385, 69328385),
  (70010542, 70060542),
  (70218729, 70318729),
  (70506535, 70556535),
  (70735468, 70835468),
  (92343416, 92443416),
  (92528796, 92678796),
  (133073060, 133223060),
  (137041193, 137091193),
  (139166997, 139216997),
  (141153431, 141203431),
  (141203431, 141213431)],
 'chrX': [(0, 10000),
  (10000, 60000),
  (94821, 144821),
  (231384, 281384),
  (1047557, 1097557),
  (1134113, 1184113),
  (1264234, 1314234),
  (2068238, 2118238),
  (7623882, 7673882),
  (10738674, 10788674),
  (37098256, 37148256),
  (49242997, 49292997),
  (49974173, 50024173),
  (52395914, 52445914),
  (58582012, 58632012),
  (58632012, 61632012),
  (61632012, 61682012),
  (76653692, 76703692),
  (113517668, 113567668),
  (115682290, 115732290),
  (120013235, 120063235),
  (143507324, 143557324),
  (148906424, 148956424),
  (149032062, 149082062),
  (152277099, 152327099),
  (155260560, 155270560)],
 'chrX_jh806590_fix': [(1587823, 1717823)],
 'chrY': [(0, 10000),
  (44821, 94821),
  (181384, 231384),
  (997557, 1047557),
  (1084113, 1134113),
  (1214234, 1264234),
  (2018238, 2068238),
  (8914955, 8964955),
  (9241322, 9291322),
  (10104553, 13104553),
  (13143954, 13193954),
  (13748578, 13798578),
  (20143885, 20193885),
  (22369679, 22419679),
  (23901428, 23951428),
  (28819361, 58819361),
  (58917656, 58967656),
  (59363566, 59373566)]}

HG38_GAPS = {'chr1': [(0, 10000),
  (177417, 227417),
  (267719, 317719),
  (471368, 521368),
  (2634220, 2684220),
  (3845268, 3995268),
  (13052998, 13102998),
  (13219912, 13319912),
  (13557162, 13607162),
  (17125658, 17175658),
  (29878082, 30028082),
  (103863906, 103913906),
  (120697156, 120747156),
  (120936695, 121086695),
  (121485434, 121535434),
  (121535434, 124535434),
  (124535434, 142535434),
  (142731022, 142781022),
  (142967761, 143117761),
  (143292816, 143342816),
  (143544525, 143644525),
  (143771002, 143871002),
  (144095783, 144145783),
  (144224481, 144274481),
  (144401744, 144451744),
  (144622413, 144672413),
  (144710724, 144810724),
  (145833118, 145883118),
  (146164650, 146214650),
  (146253299, 146303299),
  (148026038, 148176038),
  (148361358, 148511358),
  (148684147, 148734147),
  (148954460, 149004460),
  (149459645, 149509645),
  (205922707, 206072707),
  (206332221, 206482221),
  (223747846, 223797846),
  (235192211, 235242211),
  (248908210, 249058210),
  (249240621, 249250621)],
 'chr10': [(0, 10000),
  (10000, 60000),
  (17974675, 18024675),
  (38818835, 38868835),
  (39154935, 39254935),
  (39254935, 42254935),
  (42254935, 42354935),
  (42546687, 42596687),
  (46426964, 46476964),
  (47429169, 47529169),
  (47792476, 47892476),
  (48055707, 48105707),
  (49095536, 49195536),
  (51137410, 51187410),
  (51398845, 51448845),
  (125869472, 125919472),
  (128616069, 128766069),
  (133381404, 133431404),
  (133677527, 133727527),
  (135524747, 135534747)],
 'chr10_jh591181_fix': [(2010264, 2060264)],
 'chr11': [(0, 10000),
  (10000, 60000),
  (1162759, 1212759),
  (50783853, 50833853),
  (50833853, 51040853),
  (51040853, 51090853),
  (51594205, 51644205),
  (51644205, 54644205),
  (54644205, 54694205),
  (69089801, 69139801),
  (69724695, 69774695),
  (87688378, 87738378),
  (96287584, 96437584),
  (134946516, 134996516),
  (134996516, 135006516)],
 'chr11_jh806581_fix': [(499909, 549909)],
 'chr12': [(0, 10000),
  (10000, 60000),
  (95739, 145739),
  (7189876, 7239876),
  (34856694, 37856694),
  (109373470, 109423470),
  (122530623, 122580623),
  (132706992, 132806992),
  (133841895, 133851895)],
 'chr13': [(0, 10000),
  (10000, 16000000),
  (16000000, 19000000),
  (19000000, 19020000),
  (86760324, 86910324),
  (112353994, 112503994),
  (114325993, 114425993),
  (114639948, 114739948),
  (115109878, 115159878),
  (115159878, 115169878)],
 'chr14': [(0, 10000),
  (10000, 16000000),
  (16000000, 19000000),
  (107289540, 107339540),
  (107339540, 107349540)],
 'chr15': [(0, 10000),
  (10000, 17000000),
  (17000000, 20000000),
  (20894633, 20935075),
  (21398819, 21885000),
  (22212114, 22262114),
  (22596193, 22646193),
  (23514853, 23564853),
  (29159443, 29209443),
  (82829645, 82879645),
  (84984473, 85034473),
  (102521392, 102531392)],
 'chr16': [(0, 10000),
  (10000, 60000),
  (8636921, 8686921),
  (34023150, 34173150),
  (35285801, 35335801),
  (35335801, 38335801),
  (38335801, 46335801),
  (46335801, 46385801),
  (88389383, 88439383),
  (90294753, 90344753),
  (90344753, 90354753)],
 'chr17': [(296626, 396626),
  (21566608, 21666608),
  (22263006, 25263006),
  (34675848, 34725848),
  (62410760, 62460760),
  (77546461, 77596461),
  (79709049, 79759049)],
 'chr17_ctg5_hap1': [(1256794, 1306794), (1588968, 1638968)],
 'chr17_gl383561_fix': [(448214, 466467),
  (513301, 513321),
  (628875, 628895),
  (635816, 635918),
  (636913, 637464)],
 'chr17_jh720447_fix': [(190775, 240775)],
 'chr17_ke332502_fix': [(100012, 100085),
  (101792, 102078),
  (104926, 105301),
  (160048, 160456),
  (181878, 182834),
  (196662, 197232),
  (267213, 267233),
  (281523, 281724),
  (321182, 321624),
  (338198, 338308)],
 'chr18': [(0, 10000),
  (15410898, 15460898),
  (15460898, 18460898),
  (18460898, 18510898),
  (52059136, 52209136),
  (72283353, 72333353),
  (75721820, 75771820),
  (78017248, 78067248),
  (78067248, 78077248)],
 'chr18_gl383571_alt': [(44225, 94225)],
 'chr19': [(0, 10000),
  (10000, 60000),
  (7346004, 7396004),
  (8687198, 8737198),
  (20523415, 20573415),
  (24631782, 24681782),
  (24681782, 27681782),
  (27681782, 27731782),
  (59118983, 59128983)],
 'chr19_gl949746_alt': [(309793, 314793), (700652, 800652)],
 'chr19_gl949747_alt': [(195631, 256631), (394094, 399094)],
 'chr19_gl949748_alt': [(195631, 720631), (761613, 811613), (845239, 877239)],
 'chr19_gl949749_alt': [(195631, 807631), (860776, 904776)],
 'chr19_gl949750_alt': [(195631, 819631), (855325, 879325)],
 'chr19_gl949751_alt': [(195631, 705631), (789618, 815618)],
 'chr19_gl949752_alt': [(253362, 258362)],
 'chr19_gl949753_alt': [(195631, 213631), (361077, 363077), (509414, 609414)],
 'chr2': [(0, 10000),
  (3529312, 3579312),
  (5018788, 5118788),
  (16279724, 16329724),
  (21153113, 21178113),
  (87668206, 87718206),
  (89630436, 89830436),
  (90321525, 90371525),
  (90545103, 91545103),
  (91545103, 91595103),
  (92326171, 95326171),
  (110109337, 110251337),
  (149690582, 149790582),
  (234003741, 234053741),
  (239801978, 239831978),
  (240784132, 240809132),
  (243102476, 243152476),
  (243189373, 243199373)],
 'chr20': [(0, 10000),
  (10000, 60000),
  (26319569, 26369569),
  (26369569, 29369569),
  (29369569, 29419569),
  (29653908, 29803908),
  (34897085, 34947085),
  (61091437, 61141437),
  (61213369, 61263369),
  (62965520, 63015520),
  (63015520, 63025520)],
 'chr21': [(0, 10000),
  (10000, 5211193),
  (5211193, 9411193),
  (9595548, 9645548),
  (9775437, 9825437),
  (10034920, 10084920),
  (10215976, 10365976),
  (10647896, 10697896),
  (11188129, 11238129),
  (11238129, 11288129),
  (11288129, 14288129),
  (14288129, 14338129),
  (42955559, 43005559),
  (44632664, 44682664),
  (48119895, 48129895)],
 'chr22': [(0, 10000),
  (10000, 13000000),
  (13000000, 16000000),
  (16000000, 16050000),
  (16697850, 16847850),
  (20509431, 20609431),
  (50364777, 50414777),
  (51244566, 51294566),
  (51294566, 51304566)],
 'chr3': [(0, 10000),
  (10000, 60000),
  (66170270, 66270270),
  (90504854, 93504854),
  (194041961, 194047251),
  (197962430, 198012430),
  (198012430, 198022430)],
 'chr4': [(0, 10000),
  (1423146, 1478646),
  (8799203, 8818203),
  (9274642, 9324642),
  (31820917, 31837417),
  (32834638, 32840638),
  (40296396, 40297096),
  (49338941, 49488941),
  (49660117, 52660117),
  (59739333, 59789333),
  (75427379, 75452279),
  (191044276, 191144276),
  (191144276, 191154276)],
 'chr4_gl877872_fix': [(216051, 266051)],
 'chr5': [(0, 10000),
  (17530657, 17580657),
  (46405641, 49405641),
  (91636128, 91686128),
  (138787073, 138837073),
  (155138727, 155188727),
  (180905260, 180915260)],
 'chr5_gl339449_alt': [(1138692, 1238692)],
 'chr6': [(0, 10000),
  (10000, 60000),
  (58087659, 58137659),
  (58780166, 58830166),
  (58830166, 61830166),
  (61830166, 61880166),
  (62128589, 62178589),
  (95680543, 95830543),
  (157559467, 157609467),
  (157641300, 157691300),
  (167942073, 168042073),
  (170279972, 170329972),
  (171055067, 171105067),
  (171105067, 171115067)],
 'chr6_apd_hap1': [(207807, 252060),
  (327074, 356458),
  (448574, 536172),
  (588761, 692953),
  (944435, 963558),
  (970040, 1029417),
  (1173934, 1198528),
  (1307951, 1317697),
  (1341208, 1344139),
  (1451581, 1556172),
  (1582549, 1598989),
  (1736596, 1786349),
  (1833811, 1937682),
  (2009922, 2164480),
  (2169002, 2254169),
  (2359974, 2761234),
  (2821679, 2858431),
  (2894657, 3037482),
  (3122973, 3136480),
  (3180323, 3355295),
  (3370179, 3415446),
  (3491066, 3594101),
  (3638667, 3642320),
  (3696623, 3736194),
  (3762870, 3766192),
  (3790366, 3970602),
  (4220467, 4274176),
  (4390029, 4597885)],
 'chr6_dbb_hap3': [(138055, 245960),
  (622457, 640380),
  (1337861, 1344119),
  (1836690, 1837200),
  (2117295, 2119813),
  (2606257, 2723615),
  (3407000, 3481304),
  (3736386, 3772614),
  (4206542, 4226964),
  (4354126, 4376794)],
 'chr6_mann_hap4': [(145852, 184243),
  (1388812, 1438812),
  (2047420, 2056121),
  (2252213, 2277026),
  (2312291, 2324610),
  (2517052, 2540962),
  (2772006, 2810953),
  (2921096, 2929556),
  (3062107, 3151453),
  (3174271, 3227203),
  (3261361, 3374492),
  (3412258, 3450867),
  (3553455, 3557741),
  (3895092, 3926341),
  (4236950, 4245095),
  (4441658, 4480941)],
 'chr6_mcf_hap5': [(135141, 138043),
  (179429, 186757),
  (229133, 282732),
  (383988, 404415),
  (491590, 560421),
  (676599, 685588),
  (994250, 1009296),
  (1177759, 1292220),
  (1562952, 1726699),
  (1815769, 1838106),
  (2091199, 2189158),
  (2241824, 2251998),
  (2300850, 2305647),
  (2712451, 2811903),
  (2947576, 2958999),
  (3154498, 3188719),
  (3329672, 3355524),
  (3602941, 3659279),
  (3899294, 3929849),
  (4149579, 4149626),
  (4284353, 4302274),
  (4422172, 4594253)],
 'chr6_qbl_hap6': [(521900, 681210),
  (2680914, 2732064),
  (3023176, 3049607),
  (3316246, 3369039),
  (3993031, 4020006)],
 'chr6_ssto_hap7': [(861800, 941242),
  (1727717, 1744408),
  (1836632, 1838892),
  (1983561, 1997824),
  (2052735, 2087415),
  (2101474, 2113875),
  (2138491, 2303018),
  (2411247, 2461193),
  (3040184, 3054966),
  (3093568, 3154699),
  (3191196, 3214366),
  (3368931, 3445987),
  (4370564, 4420564),
  (4580410, 4610218),
  (4639806, 4661556),
  (4723509, 4774367),
  (4783798, 4836049)],
 'chr7': [(0, 10000),
  (232484, 282484),
  (50370631, 50410631),
  (58054331, 61054331),
  (61310513, 61360513),
  (61460465, 61510465),
  (61677020, 61727020),
  (61917157, 61967157),
  (74715724, 74765724),
  (100556043, 100606043),
  (130154523, 130254523),
  (139379377, 139404377),
  (142048195, 142098195),
  (142276197, 142326197),
  (143347897, 143397897),
  (154270634, 154370634),
  (159128663, 159138663)],
 'chr7_jh636058_fix': [(170100, 220100)],
 'chr8': [(0, 10000),
  (7474649, 7524649),
  (12091854, 12141854),
  (43838887, 46838887),
  (48130499, 48135599),
  (86576451, 86726451),
  (142766515, 142816515),
  (145332588, 145432588),
  (146304022, 146354022),
  (146354022, 146364022)],
 'chr9': [(0, 10000),
  (39663686, 39713686),
  (39974796, 40024796),
  (40233029, 40283029),
  (40425834, 40475834),
  (40940341, 40990341),
  (41143214, 41193214),
  (41365793, 41415793),
  (42613955, 42663955),
  (43213698, 43313698),
  (43946569, 43996569),
  (44676646, 44726646),
  (44908293, 44958293),
  (45250203, 45350203),
  (45815521, 45865521),
  (46216430, 46266430),
  (46461039, 46561039),
  (47060133, 47160133),
  (47317679, 47367679),
  (47367679, 50367679),
  (50367679, 65367679),
  (65367679, 65467679),
  (65918360, 65968360),
  (66192215, 66242215),
  (66404656, 66454656),
  (66614195, 66664195),
  (66863343, 66913343),
  (67107834, 67207834),
  (67366296, 67516296),
  (67987998, 68137998),
  (68514181, 68664181),
  (68838946, 68988946),
  (69278385, 69328385),
  (70010542, 70060542),
  (70218729, 70318729),
  (70506535, 70556535),
  (70735468, 70835468),
  (92343416, 92443416),
  (92528796, 92678796),
  (133073060, 133223060),
  (137041193, 137091193),
  (139166997, 139216997),
  (141153431, 141203431),
  (141203431, 141213431)],
 'chrX': [(0, 10000),
  (10000, 60000),
  (94821, 144821),
  (231384, 281384),
  (1047557, 1097557),
  (1134113, 1184113),
  (1264234, 1314234),
  (2068238, 2118238),
  (7623882, 7673882),
  (10738674, 10788674),
  (37098256, 37148256),
  (49242997, 49292997),
  (49974173, 50024173),
  (52395914, 52445914),
  (58582012, 58632012),
  (58632012, 61632012),
  (61632012, 61682012),
  (76653692, 76703692),
  (113517668, 113567668),
  (115682290, 115732290),
  (120013235, 120063235),
  (143507324, 143557324),
  (148906424, 148956424),
  (149032062, 149082062),
  (152277099, 152327099),
  (155260560, 155270560)],
 'chrX_jh806590_fix': [(1587823, 1717823)],
 'chrY': [(0, 10000),
  (44821, 94821),
  (181384, 231384),
  (997557, 1047557),
  (1084113, 1134113),
  (1214234, 1264234),
  (2018238, 2068238),
  (8914955, 8964955),
  (9241322, 9291322),
  (10104553, 13104553),
  (13143954, 13193954),
  (13748578, 13798578),
  (20143885, 20193885),
  (22369679, 22419679),
  (23901428, 23951428),
  (28819361, 58819361),
  (58917656, 58967656),
  (59363566, 59373566)]}
