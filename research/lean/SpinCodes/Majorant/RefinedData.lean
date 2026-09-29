import SpinCodes.Structured.Majorant
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace Spin.Majorant
open Spin.Structured
def supports : Finset (ℚ × ℚ) := {((833 / 500), (-43311 / 250000)), ((83 / 50), (-3446583973879 / 20000000000000)), ((33 / 20), (-5335812508647 / 31250000000000)), ((163 / 100), (-3350558688107 / 20000000000000)), ((8 / 5), (-40647089344673 / 250000000000000)), ((31 / 20), (-77022601619127 / 500000000000000)), ((7 / 5), (-31533292681423 / 250000000000000)), ((6 / 5), (-832718284462267 / 10000000000000000)), ((1 / 1), (-16647078810249 / 500000000000000)), ((4 / 5), (30683842683431 / 1250000000000000)), ((3 / 5), (909341020092481 / 10000000000000000)), ((2 / 5), (4161516955371 / 25000000000000)), ((1 / 5), (15724076464601 / 62500000000000)), ((1 / 10), (297841147503783 / 1000000000000000)), ((9 / 100), (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456)), ((7 / 100), (53117391002404332784602253741220096935 / 170141183460469231731687303715884105728)), ((1 / 20), (109535552428011023053662565657799459079 / 340282366920938463463374607431768211456)), ((3 / 100), (112870329263410277722384480428476839213 / 340282366920938463463374607431768211456)), ((1 / 100), (58119563056842377037993588567825857659 / 170141183460469231731687303715884105728)), ((0 / 1), (3465735902799727 / 10000000000000000)), ((-1 / 100), (1495524372286176733882761540124617467907 / 4253529586511730793292182592897102643200)), ((-3 / 100), (3076970006775960790657142966285747138917 / 8507059173023461586584365185794205286400)), ((-1 / 20), (632748353870289731134156480146939348259 / 1701411834604692317316873037158841057280)), ((-7 / 100), (1625681846115929475145509125033299608399 / 4253529586511730793292182592897102643200)), ((-9 / 100), (3339835265656588966103687050071259581601 / 8507059173023461586584365185794205286400)), ((-1 / 10), (397841147503783 / 1000000000000000)), ((-1 / 5), (28224076464601 / 62500000000000)), ((-2 / 5), (14161516955371 / 25000000000000)), ((-3 / 5), (6909341020092481 / 10000000000000000)), ((-4 / 5), (1030683842683431 / 1250000000000000)), ((-1 / 1), (483352921189751 / 500000000000000)), ((-6 / 5), (11167281715537733 / 10000000000000000)), ((-7 / 5), (318466707318577 / 250000000000000)), ((-31 / 20), (697977398380873 / 500000000000000)), ((-8 / 5), (359352910655327 / 250000000000000)), ((-163 / 100), (29249441311893 / 20000000000000)), ((-33 / 20), (46226687491353 / 31250000000000)), ((-83 / 50), (29753416026121 / 20000000000000)), ((-833 / 500), (373189 / 250000))}
def refined : ConcaveMajorant := ⟨supports, by exact ⟨_, Finset.mem_insert_self _ _⟩⟩
def reflect (p : ℚ × ℚ) : ℚ × ℚ := (-p.1, p.1+p.2)
lemma member0 : ((833 / 500), (-43311 / 250000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inl True.intro)
lemma member1 : ((83 / 50), (-3446583973879 / 20000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inl True.intro))
lemma member2 : ((33 / 20), (-5335812508647 / 31250000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inl True.intro)))
lemma member3 : ((163 / 100), (-3350558688107 / 20000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))
lemma member4 : ((8 / 5), (-40647089344673 / 250000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))
lemma member5 : ((31 / 20), (-77022601619127 / 500000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))
lemma member6 : ((7 / 5), (-31533292681423 / 250000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))
lemma member7 : ((6 / 5), (-832718284462267 / 10000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))
lemma member8 : ((1 / 1), (-16647078810249 / 500000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))
lemma member9 : ((4 / 5), (30683842683431 / 1250000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))
lemma member10 : ((3 / 5), (909341020092481 / 10000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))
lemma member11 : ((2 / 5), (4161516955371 / 25000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))
lemma member12 : ((1 / 5), (15724076464601 / 62500000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))
lemma member13 : ((1 / 10), (297841147503783 / 1000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))
lemma member14 : ((9 / 100), (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))
lemma member15 : ((7 / 100), (53117391002404332784602253741220096935 / 170141183460469231731687303715884105728)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))
lemma member16 : ((1 / 20), (109535552428011023053662565657799459079 / 340282366920938463463374607431768211456)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))
lemma member17 : ((3 / 100), (112870329263410277722384480428476839213 / 340282366920938463463374607431768211456)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))
lemma member18 : ((1 / 100), (58119563056842377037993588567825857659 / 170141183460469231731687303715884105728)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))
lemma member19 : ((0 / 1), (3465735902799727 / 10000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))
lemma member20 : ((-1 / 100), (1495524372286176733882761540124617467907 / 4253529586511730793292182592897102643200)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))
lemma member21 : ((-3 / 100), (3076970006775960790657142966285747138917 / 8507059173023461586584365185794205286400)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))
lemma member22 : ((-1 / 20), (632748353870289731134156480146939348259 / 1701411834604692317316873037158841057280)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))
lemma member23 : ((-7 / 100), (1625681846115929475145509125033299608399 / 4253529586511730793292182592897102643200)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))
lemma member24 : ((-9 / 100), (3339835265656588966103687050071259581601 / 8507059173023461586584365185794205286400)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))
lemma member25 : ((-1 / 10), (397841147503783 / 1000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))
lemma member26 : ((-1 / 5), (28224076464601 / 62500000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))))
lemma member27 : ((-2 / 5), (14161516955371 / 25000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))))
lemma member28 : ((-3 / 5), (6909341020092481 / 10000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))))))
lemma member29 : ((-4 / 5), (1030683842683431 / 1250000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))))))
lemma member30 : ((-1 / 1), (483352921189751 / 500000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))))))))
lemma member31 : ((-6 / 5), (11167281715537733 / 10000000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))))))))
lemma member32 : ((-7 / 5), (318466707318577 / 250000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))))))))))
lemma member33 : ((-31 / 20), (697977398380873 / 500000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))))))))))
lemma member34 : ((-8 / 5), (359352910655327 / 250000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))))))))))))
lemma member35 : ((-163 / 100), (29249441311893 / 20000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))))))))))))
lemma member36 : ((-33 / 20), (46226687491353 / 31250000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro)))))))))))))))))))))))))))))))))))))
lemma member37 : ((-83 / 50), (29753416026121 / 20000000000000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl True.intro))))))))))))))))))))))))))))))))))))))
lemma member38 : ((-833 / 500), (373189 / 250000)) ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr True.intro))))))))))))))))))))))))))))))))))))))
theorem reflect_mem {p : ℚ × ℚ} (hp : p ∈ supports) : reflect p ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton] at hp
  rcases hp with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · have he : reflect ((833 / 500), (-43311 / 250000)) = ((-833 / 500), (373189 / 250000)) := by norm_num [reflect]
    rw [he]
    exact member38
  · have he : reflect ((83 / 50), (-3446583973879 / 20000000000000)) = ((-83 / 50), (29753416026121 / 20000000000000)) := by norm_num [reflect]
    rw [he]
    exact member37
  · have he : reflect ((33 / 20), (-5335812508647 / 31250000000000)) = ((-33 / 20), (46226687491353 / 31250000000000)) := by norm_num [reflect]
    rw [he]
    exact member36
  · have he : reflect ((163 / 100), (-3350558688107 / 20000000000000)) = ((-163 / 100), (29249441311893 / 20000000000000)) := by norm_num [reflect]
    rw [he]
    exact member35
  · have he : reflect ((8 / 5), (-40647089344673 / 250000000000000)) = ((-8 / 5), (359352910655327 / 250000000000000)) := by norm_num [reflect]
    rw [he]
    exact member34
  · have he : reflect ((31 / 20), (-77022601619127 / 500000000000000)) = ((-31 / 20), (697977398380873 / 500000000000000)) := by norm_num [reflect]
    rw [he]
    exact member33
  · have he : reflect ((7 / 5), (-31533292681423 / 250000000000000)) = ((-7 / 5), (318466707318577 / 250000000000000)) := by norm_num [reflect]
    rw [he]
    exact member32
  · have he : reflect ((6 / 5), (-832718284462267 / 10000000000000000)) = ((-6 / 5), (11167281715537733 / 10000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member31
  · have he : reflect ((1 / 1), (-16647078810249 / 500000000000000)) = ((-1 / 1), (483352921189751 / 500000000000000)) := by norm_num [reflect]
    rw [he]
    exact member30
  · have he : reflect ((4 / 5), (30683842683431 / 1250000000000000)) = ((-4 / 5), (1030683842683431 / 1250000000000000)) := by norm_num [reflect]
    rw [he]
    exact member29
  · have he : reflect ((3 / 5), (909341020092481 / 10000000000000000)) = ((-3 / 5), (6909341020092481 / 10000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member28
  · have he : reflect ((2 / 5), (4161516955371 / 25000000000000)) = ((-2 / 5), (14161516955371 / 25000000000000)) := by norm_num [reflect]
    rw [he]
    exact member27
  · have he : reflect ((1 / 5), (15724076464601 / 62500000000000)) = ((-1 / 5), (28224076464601 / 62500000000000)) := by norm_num [reflect]
    rw [he]
    exact member26
  · have he : reflect ((1 / 10), (297841147503783 / 1000000000000000)) = ((-1 / 10), (397841147503783 / 1000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member25
  · have he : reflect ((9 / 100), (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456)) = ((-9 / 100), (3339835265656588966103687050071259581601 / 8507059173023461586584365185794205286400)) := by norm_num [reflect]
    rw [he]
    exact member24
  · have he : reflect ((7 / 100), (53117391002404332784602253741220096935 / 170141183460469231731687303715884105728)) = ((-7 / 100), (1625681846115929475145509125033299608399 / 4253529586511730793292182592897102643200)) := by norm_num [reflect]
    rw [he]
    exact member23
  · have he : reflect ((1 / 20), (109535552428011023053662565657799459079 / 340282366920938463463374607431768211456)) = ((-1 / 20), (632748353870289731134156480146939348259 / 1701411834604692317316873037158841057280)) := by norm_num [reflect]
    rw [he]
    exact member22
  · have he : reflect ((3 / 100), (112870329263410277722384480428476839213 / 340282366920938463463374607431768211456)) = ((-3 / 100), (3076970006775960790657142966285747138917 / 8507059173023461586584365185794205286400)) := by norm_num [reflect]
    rw [he]
    exact member21
  · have he : reflect ((1 / 100), (58119563056842377037993588567825857659 / 170141183460469231731687303715884105728)) = ((-1 / 100), (1495524372286176733882761540124617467907 / 4253529586511730793292182592897102643200)) := by norm_num [reflect]
    rw [he]
    exact member20
  · have he : reflect ((0 / 1), (3465735902799727 / 10000000000000000)) = ((0 / 1), (3465735902799727 / 10000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member19
  · have he : reflect ((-1 / 100), (1495524372286176733882761540124617467907 / 4253529586511730793292182592897102643200)) = ((1 / 100), (58119563056842377037993588567825857659 / 170141183460469231731687303715884105728)) := by norm_num [reflect]
    rw [he]
    exact member18
  · have he : reflect ((-3 / 100), (3076970006775960790657142966285747138917 / 8507059173023461586584365185794205286400)) = ((3 / 100), (112870329263410277722384480428476839213 / 340282366920938463463374607431768211456)) := by norm_num [reflect]
    rw [he]
    exact member17
  · have he : reflect ((-1 / 20), (632748353870289731134156480146939348259 / 1701411834604692317316873037158841057280)) = ((1 / 20), (109535552428011023053662565657799459079 / 340282366920938463463374607431768211456)) := by norm_num [reflect]
    rw [he]
    exact member16
  · have he : reflect ((-7 / 100), (1625681846115929475145509125033299608399 / 4253529586511730793292182592897102643200)) = ((7 / 100), (53117391002404332784602253741220096935 / 170141183460469231731687303715884105728)) := by norm_num [reflect]
    rw [he]
    exact member15
  · have he : reflect ((-9 / 100), (3339835265656588966103687050071259581601 / 8507059173023461586584365185794205286400)) = ((9 / 100), (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456)) := by norm_num [reflect]
    rw [he]
    exact member14
  · have he : reflect ((-1 / 10), (397841147503783 / 1000000000000000)) = ((1 / 10), (297841147503783 / 1000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member13
  · have he : reflect ((-1 / 5), (28224076464601 / 62500000000000)) = ((1 / 5), (15724076464601 / 62500000000000)) := by norm_num [reflect]
    rw [he]
    exact member12
  · have he : reflect ((-2 / 5), (14161516955371 / 25000000000000)) = ((2 / 5), (4161516955371 / 25000000000000)) := by norm_num [reflect]
    rw [he]
    exact member11
  · have he : reflect ((-3 / 5), (6909341020092481 / 10000000000000000)) = ((3 / 5), (909341020092481 / 10000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member10
  · have he : reflect ((-4 / 5), (1030683842683431 / 1250000000000000)) = ((4 / 5), (30683842683431 / 1250000000000000)) := by norm_num [reflect]
    rw [he]
    exact member9
  · have he : reflect ((-1 / 1), (483352921189751 / 500000000000000)) = ((1 / 1), (-16647078810249 / 500000000000000)) := by norm_num [reflect]
    rw [he]
    exact member8
  · have he : reflect ((-6 / 5), (11167281715537733 / 10000000000000000)) = ((6 / 5), (-832718284462267 / 10000000000000000)) := by norm_num [reflect]
    rw [he]
    exact member7
  · have he : reflect ((-7 / 5), (318466707318577 / 250000000000000)) = ((7 / 5), (-31533292681423 / 250000000000000)) := by norm_num [reflect]
    rw [he]
    exact member6
  · have he : reflect ((-31 / 20), (697977398380873 / 500000000000000)) = ((31 / 20), (-77022601619127 / 500000000000000)) := by norm_num [reflect]
    rw [he]
    exact member5
  · have he : reflect ((-8 / 5), (359352910655327 / 250000000000000)) = ((8 / 5), (-40647089344673 / 250000000000000)) := by norm_num [reflect]
    rw [he]
    exact member4
  · have he : reflect ((-163 / 100), (29249441311893 / 20000000000000)) = ((163 / 100), (-3350558688107 / 20000000000000)) := by norm_num [reflect]
    rw [he]
    exact member3
  · have he : reflect ((-33 / 20), (46226687491353 / 31250000000000)) = ((33 / 20), (-5335812508647 / 31250000000000)) := by norm_num [reflect]
    rw [he]
    exact member2
  · have he : reflect ((-83 / 50), (29753416026121 / 20000000000000)) = ((83 / 50), (-3446583973879 / 20000000000000)) := by norm_num [reflect]
    rw [he]
    exact member1
  · have he : reflect ((-833 / 500), (373189 / 250000)) = ((833 / 500), (-43311 / 250000)) := by norm_num [reflect]
    rw [he]
    exact member0

theorem refined_reflect (w : ℝ) : refined.toFun (1-w) = refined.toFun w := by
  apply le_antisymm
  · apply Finset.le_inf'
    intro p hp
    have h := refined.toFun_le_support (reflect_mem hp) (1-w)
    convert h using 1 <;> simp only [reflect, Rat.cast_neg, Rat.cast_add] <;> ring
  · apply Finset.le_inf'
    intro p hp
    have h := refined.toFun_le_support (reflect_mem hp) w
    convert h using 1 <;> simp only [reflect, Rat.cast_neg, Rat.cast_add] <;> ring

lemma affine_le_between {s i t j l u x : ℝ} (hx : x ∈ Set.Icc l u)
    (hl : s*l+i ≤ t*l+j) (hu : s*u+i ≤ t*u+j) : s*x+i ≤ t*x+j := by
  rcases le_total s t with h | h
  · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hx.1)]
  · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hx.2)]

lemma active_le_refined {s i l u : ℚ} {w : ℝ} (hw : w ∈ Set.Icc (l : ℝ) (u : ℝ))
    (h : ∀ p ∈ supports, s*l+i ≤ p.1*l+p.2 ∧ s*u+i ≤ p.1*u+p.2) :
    (s : ℝ)*w+(i : ℝ) ≤ refined.toFun w := by
  apply Finset.le_inf'
  intro p hp
  obtain ⟨hl, hu⟩ := h p hp
  apply affine_le_between hw
  · exact_mod_cast hl
  · exact_mod_cast hu
theorem active0 : ∀ p ∈ supports,
    ((833 / 500):ℚ)*(13 / 125)+(-43311 / 250000) ≤ p.1*(13 / 125)+p.2 ∧
    ((833 / 500):ℚ)*(18296026121 / 120000000000)+(-43311 / 250000) ≤ p.1*(18296026121 / 120000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active1 : ∀ p ∈ supports,
    ((83 / 50):ℚ)*(18296026121 / 120000000000)+(-3446583973879 / 20000000000000) ≤ p.1*(18296026121 / 120000000000)+p.2 ∧
    ((83 / 50):ℚ)*(791599208623 / 5000000000000)+(-3446583973879 / 20000000000000) ≤ p.1*(791599208623 / 5000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active2 : ∀ p ∈ supports,
    ((33 / 20):ℚ)*(791599208623 / 5000000000000)+(-5335812508647 / 31250000000000) ≤ p.1*(791599208623 / 5000000000000)+p.2 ∧
    ((33 / 20):ℚ)*(1609032935677 / 10000000000000)+(-5335812508647 / 31250000000000) ≤ p.1*(1609032935677 / 10000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active3 : ∀ p ∈ supports,
    ((163 / 100):ℚ)*(1609032935677 / 10000000000000)+(-3350558688107 / 20000000000000) ≤ p.1*(1609032935677 / 10000000000000)+p.2 ∧
    ((163 / 100):ℚ)*(2469788513329 / 15000000000000)+(-3350558688107 / 20000000000000) ≤ p.1*(2469788513329 / 15000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active4 : ∀ p ∈ supports,
    ((8 / 5):ℚ)*(2469788513329 / 15000000000000)+(-40647089344673 / 250000000000000) ≤ p.1*(2469788513329 / 15000000000000)+p.2 ∧
    ((8 / 5):ℚ)*(4271577070219 / 25000000000000)+(-40647089344673 / 250000000000000) ≤ p.1*(4271577070219 / 25000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active5 : ∀ p ∈ supports,
    ((31 / 20):ℚ)*(4271577070219 / 25000000000000)+(-77022601619127 / 500000000000000) ≤ p.1*(4271577070219 / 25000000000000)+p.2 ∧
    ((31 / 20):ℚ)*(13956016256281 / 75000000000000)+(-77022601619127 / 500000000000000) ≤ p.1*(13956016256281 / 75000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active6 : ∀ p ∈ supports,
    ((7 / 5):ℚ)*(13956016256281 / 75000000000000)+(-31533292681423 / 250000000000000) ≤ p.1*(13956016256281 / 75000000000000)+p.2 ∧
    ((7 / 5):ℚ)*(428613422794653 / 2000000000000000)+(-31533292681423 / 250000000000000) ≤ p.1*(428613422794653 / 2000000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active7 : ∀ p ∈ supports,
    ((6 / 5):ℚ)*(428613422794653 / 2000000000000000)+(-832718284462267 / 10000000000000000) ≤ p.1*(428613422794653 / 2000000000000000)+p.2 ∧
    ((6 / 5):ℚ)*(499776708257287 / 2000000000000000)+(-832718284462267 / 10000000000000000) ≤ p.1*(499776708257287 / 2000000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active8 : ∀ p ∈ supports,
    ((1 / 1):ℚ)*(499776708257287 / 2000000000000000)+(-16647078810249 / 500000000000000) ≤ p.1*(499776708257287 / 2000000000000000)+p.2 ∧
    ((1 / 1):ℚ)*(144603079418107 / 500000000000000)+(-16647078810249 / 500000000000000) ≤ p.1*(144603079418107 / 500000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active9 : ∀ p ∈ supports,
    ((4 / 5):ℚ)*(144603079418107 / 500000000000000)+(30683842683431 / 1250000000000000) ≤ p.1*(144603079418107 / 500000000000000)+p.2 ∧
    ((4 / 5):ℚ)*(663870278625033 / 2000000000000000)+(30683842683431 / 1250000000000000) ≤ p.1*(663870278625033 / 2000000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active10 : ∀ p ∈ supports,
    ((3 / 5):ℚ)*(663870278625033 / 2000000000000000)+(909341020092481 / 10000000000000000) ≤ p.1*(663870278625033 / 2000000000000000)+p.2 ∧
    ((3 / 5):ℚ)*(755265762055919 / 2000000000000000)+(909341020092481 / 10000000000000000) ≤ p.1*(755265762055919 / 2000000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active11 : ∀ p ∈ supports,
    ((2 / 5):ℚ)*(755265762055919 / 2000000000000000)+(4161516955371 / 25000000000000) ≤ p.1*(755265762055919 / 2000000000000000)+p.2 ∧
    ((2 / 5):ℚ)*(10640568152347 / 25000000000000)+(4161516955371 / 25000000000000) ≤ p.1*(10640568152347 / 25000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active12 : ∀ p ∈ supports,
    ((1 / 5):ℚ)*(10640568152347 / 25000000000000)+(15724076464601 / 62500000000000) ≤ p.1*(10640568152347 / 25000000000000)+p.2 ∧
    ((1 / 5):ℚ)*(46255924070167 / 100000000000000)+(15724076464601 / 62500000000000) ≤ p.1*(46255924070167 / 100000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active13 : ∀ p ∈ supports,
    ((1 / 10):ℚ)*(46255924070167 / 100000000000000)+(297841147503783 / 1000000000000000) ≤ p.1*(46255924070167 / 100000000000000)+p.2 ∧
    ((1 / 10):ℚ)*(49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000)+(297841147503783 / 1000000000000000) ≤ p.1*(49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active14 : ∀ p ∈ supports,
    ((9 / 100):ℚ)*(49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000)+(102967997603379096932443767333991244233 / 340282366920938463463374607431768211456) ≤ p.1*(49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000)+p.2 ∧
    ((9 / 100):ℚ)*(81669610035739215919018503711223740925 / 170141183460469231731687303715884105728)+(102967997603379096932443767333991244233 / 340282366920938463463374607431768211456) ≤ p.1*(81669610035739215919018503711223740925 / 170141183460469231731687303715884105728)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active15 : ∀ p ∈ supports,
    ((7 / 100):ℚ)*(81669610035739215919018503711223740925 / 170141183460469231731687303715884105728)+(53117391002404332784602253741220096935 / 170141183460469231731687303715884105728) ≤ p.1*(81669610035739215919018503711223740925 / 170141183460469231731687303715884105728)+p.2 ∧
    ((7 / 100):ℚ)*(82519260580058937111451454383981630225 / 170141183460469231731687303715884105728)+(53117391002404332784602253741220096935 / 170141183460469231731687303715884105728) ≤ p.1*(82519260580058937111451454383981630225 / 170141183460469231731687303715884105728)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active16 : ∀ p ∈ supports,
    ((1 / 20):ℚ)*(82519260580058937111451454383981630225 / 170141183460469231731687303715884105728)+(109535552428011023053662565657799459079 / 340282366920938463463374607431768211456) ≤ p.1*(82519260580058937111451454383981630225 / 170141183460469231731687303715884105728)+p.2 ∧
    ((1 / 20):ℚ)*(41684710442490683359023934633467251675 / 85070591730234615865843651857942052864)+(109535552428011023053662565657799459079 / 340282366920938463463374607431768211456) ≤ p.1*(41684710442490683359023934633467251675 / 85070591730234615865843651857942052864)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active17 : ∀ p ∈ supports,
    ((3 / 100):ℚ)*(41684710442490683359023934633467251675 / 85070591730234615865843651857942052864)+(112870329263410277722384480428476839213 / 340282366920938463463374607431768211456) ≤ p.1*(41684710442490683359023934633467251675 / 85070591730234615865843651857942052864)+p.2 ∧
    ((3 / 100):ℚ)*(84219921256861908840067417679371902625 / 170141183460469231731687303715884105728)+(112870329263410277722384480428476839213 / 340282366920938463463374607431768211456) ≤ p.1*(84219921256861908840067417679371902625 / 170141183460469231731687303715884105728)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active18 : ∀ p ∈ supports,
    ((1 / 100):ℚ)*(84219921256861908840067417679371902625 / 170141183460469231731687303715884105728)+(58119563056842377037993588567825857659 / 170141183460469231731687303715884105728) ≤ p.1*(84219921256861908840067417679371902625 / 170141183460469231731687303715884105728)+p.2 ∧
    ((1 / 100):ℚ)*(129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000)+(58119563056842377037993588567825857659 / 170141183460469231731687303715884105728) ≤ p.1*(129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active19 : ∀ p ∈ supports,
    ((0 / 1):ℚ)*(129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000)+(3465735902799727 / 10000000000000000) ≤ p.1*(129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000)+p.2 ∧
    ((0 / 1):ℚ)*(1 / 2)+(3465735902799727 / 10000000000000000) ≤ p.1*(1 / 2)+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
theorem active_real0 {w : ℝ} (hw : w ∈ Set.Icc ((13 / 125):ℝ) (18296026121 / 120000000000)) :
    ((833 / 500):ℝ)*w+(-43311 / 250000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((13 / 125):ℚ):ℝ) (((18296026121 / 120000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active0
  norm_num at h
  linarith
theorem active_real1 {w : ℝ} (hw : w ∈ Set.Icc ((18296026121 / 120000000000):ℝ) (791599208623 / 5000000000000)) :
    ((83 / 50):ℝ)*w+(-3446583973879 / 20000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((18296026121 / 120000000000):ℚ):ℝ) (((791599208623 / 5000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active1
  norm_num at h
  linarith
theorem active_real2 {w : ℝ} (hw : w ∈ Set.Icc ((791599208623 / 5000000000000):ℝ) (1609032935677 / 10000000000000)) :
    ((33 / 20):ℝ)*w+(-5335812508647 / 31250000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((791599208623 / 5000000000000):ℚ):ℝ) (((1609032935677 / 10000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active2
  norm_num at h
  linarith
theorem active_real3 {w : ℝ} (hw : w ∈ Set.Icc ((1609032935677 / 10000000000000):ℝ) (2469788513329 / 15000000000000)) :
    ((163 / 100):ℝ)*w+(-3350558688107 / 20000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((1609032935677 / 10000000000000):ℚ):ℝ) (((2469788513329 / 15000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active3
  norm_num at h
  linarith
theorem active_real4 {w : ℝ} (hw : w ∈ Set.Icc ((2469788513329 / 15000000000000):ℝ) (4271577070219 / 25000000000000)) :
    ((8 / 5):ℝ)*w+(-40647089344673 / 250000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((2469788513329 / 15000000000000):ℚ):ℝ) (((4271577070219 / 25000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active4
  norm_num at h
  linarith
theorem active_real5 {w : ℝ} (hw : w ∈ Set.Icc ((4271577070219 / 25000000000000):ℝ) (13956016256281 / 75000000000000)) :
    ((31 / 20):ℝ)*w+(-77022601619127 / 500000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((4271577070219 / 25000000000000):ℚ):ℝ) (((13956016256281 / 75000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active5
  norm_num at h
  linarith
theorem active_real6 {w : ℝ} (hw : w ∈ Set.Icc ((13956016256281 / 75000000000000):ℝ) (428613422794653 / 2000000000000000)) :
    ((7 / 5):ℝ)*w+(-31533292681423 / 250000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((13956016256281 / 75000000000000):ℚ):ℝ) (((428613422794653 / 2000000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active6
  norm_num at h
  linarith
theorem active_real7 {w : ℝ} (hw : w ∈ Set.Icc ((428613422794653 / 2000000000000000):ℝ) (499776708257287 / 2000000000000000)) :
    ((6 / 5):ℝ)*w+(-832718284462267 / 10000000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((428613422794653 / 2000000000000000):ℚ):ℝ) (((499776708257287 / 2000000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active7
  norm_num at h
  linarith
theorem active_real8 {w : ℝ} (hw : w ∈ Set.Icc ((499776708257287 / 2000000000000000):ℝ) (144603079418107 / 500000000000000)) :
    ((1 / 1):ℝ)*w+(-16647078810249 / 500000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((499776708257287 / 2000000000000000):ℚ):ℝ) (((144603079418107 / 500000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active8
  norm_num at h
  linarith
theorem active_real9 {w : ℝ} (hw : w ∈ Set.Icc ((144603079418107 / 500000000000000):ℝ) (663870278625033 / 2000000000000000)) :
    ((4 / 5):ℝ)*w+(30683842683431 / 1250000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((144603079418107 / 500000000000000):ℚ):ℝ) (((663870278625033 / 2000000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active9
  norm_num at h
  linarith
theorem active_real10 {w : ℝ} (hw : w ∈ Set.Icc ((663870278625033 / 2000000000000000):ℝ) (755265762055919 / 2000000000000000)) :
    ((3 / 5):ℝ)*w+(909341020092481 / 10000000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((663870278625033 / 2000000000000000):ℚ):ℝ) (((755265762055919 / 2000000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active10
  norm_num at h
  linarith
theorem active_real11 {w : ℝ} (hw : w ∈ Set.Icc ((755265762055919 / 2000000000000000):ℝ) (10640568152347 / 25000000000000)) :
    ((2 / 5):ℝ)*w+(4161516955371 / 25000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((755265762055919 / 2000000000000000):ℚ):ℝ) (((10640568152347 / 25000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active11
  norm_num at h
  linarith
theorem active_real12 {w : ℝ} (hw : w ∈ Set.Icc ((10640568152347 / 25000000000000):ℝ) (46255924070167 / 100000000000000)) :
    ((1 / 5):ℝ)*w+(15724076464601 / 62500000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((10640568152347 / 25000000000000):ℚ):ℝ) (((46255924070167 / 100000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active12
  norm_num at h
  linarith
theorem active_real13 {w : ℝ} (hw : w ∈ Set.Icc ((46255924070167 / 100000000000000):ℝ) (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000)) :
    ((1 / 10):ℝ)*w+(297841147503783 / 1000000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((46255924070167 / 100000000000000):ℚ):ℝ) (((49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active13
  norm_num at h
  linarith
theorem active_real14 {w : ℝ} (hw : w ∈ Set.Icc ((49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000):ℝ) (81669610035739215919018503711223740925 / 170141183460469231731687303715884105728)) :
    ((9 / 100):ℝ)*w+(102967997603379096932443767333991244233 / 340282366920938463463374607431768211456) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000):ℚ):ℝ) (((81669610035739215919018503711223740925 / 170141183460469231731687303715884105728):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active14
  norm_num at h
  linarith
theorem active_real15 {w : ℝ} (hw : w ∈ Set.Icc ((81669610035739215919018503711223740925 / 170141183460469231731687303715884105728):ℝ) (82519260580058937111451454383981630225 / 170141183460469231731687303715884105728)) :
    ((7 / 100):ℝ)*w+(53117391002404332784602253741220096935 / 170141183460469231731687303715884105728) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((81669610035739215919018503711223740925 / 170141183460469231731687303715884105728):ℚ):ℝ) (((82519260580058937111451454383981630225 / 170141183460469231731687303715884105728):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active15
  norm_num at h
  linarith
theorem active_real16 {w : ℝ} (hw : w ∈ Set.Icc ((82519260580058937111451454383981630225 / 170141183460469231731687303715884105728):ℝ) (41684710442490683359023934633467251675 / 85070591730234615865843651857942052864)) :
    ((1 / 20):ℝ)*w+(109535552428011023053662565657799459079 / 340282366920938463463374607431768211456) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((82519260580058937111451454383981630225 / 170141183460469231731687303715884105728):ℚ):ℝ) (((41684710442490683359023934633467251675 / 85070591730234615865843651857942052864):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active16
  norm_num at h
  linarith
theorem active_real17 {w : ℝ} (hw : w ∈ Set.Icc ((41684710442490683359023934633467251675 / 85070591730234615865843651857942052864):ℝ) (84219921256861908840067417679371902625 / 170141183460469231731687303715884105728)) :
    ((3 / 100):ℝ)*w+(112870329263410277722384480428476839213 / 340282366920938463463374607431768211456) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((41684710442490683359023934633467251675 / 85070591730234615865843651857942052864):ℚ):ℝ) (((84219921256861908840067417679371902625 / 170141183460469231731687303715884105728):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active17
  norm_num at h
  linarith
theorem active_real18 {w : ℝ} (hw : w ∈ Set.Icc ((84219921256861908840067417679371902625 / 170141183460469231731687303715884105728):ℝ) (129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000)) :
    ((1 / 100):ℝ)*w+(58119563056842377037993588567825857659 / 170141183460469231731687303715884105728) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((84219921256861908840067417679371902625 / 170141183460469231731687303715884105728):ℚ):ℝ) (((129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active18
  norm_num at h
  linarith
theorem active_real19 {w : ℝ} (hw : w ∈ Set.Icc ((129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000):ℝ) (1 / 2)) :
    ((0 / 1):ℝ)*w+(3465735902799727 / 10000000000000000) ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (((129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000):ℚ):ℝ) (((1 / 2):ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active19
  norm_num at h
  linarith
end Spin.Majorant
