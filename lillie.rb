# before first run:  $ bundle
# run             :  $ bundle exec ruby lillie.rb


# 2016-12-15 @yuaaa

require 'bundler'
require 'minisat'
@solver = MiniSat::Solver.new

MAX_SUTE1 = 4
MAX_SUTE2 = 4


# > 先親遺伝
# < 後親遺伝
# _ 破棄（一個目決定後）
# | 破棄（二個目決定後）
# - 偶数個体値
# + 奇数個体値

## 偶然の一致でランダム生成された個体値が片親と同じになるのは構わないが、
## 両親の個体値はそれぞれ違うものにしておくこと

oya  = []
oya[0]  = [23, 13, 18, 13, 20,  6]
oya[1]  = [30, 17, 26,  8, 29 ,13]

poke = []
poke[0] = [ 8, 17,  4, 13, 23, 13]
poke[1] = [19, 17, 31, 23, 20, 13]
poke[2] = [30, 31, 23, 13,  7, 13]
poke[3] = [30, 17, 23, 13,  7,  6]
poke[4] = [31, 23,  2, 13, 20, 13]
poke[5] = [30,  2, 26, 13, 19, 10]

def pretty(str)
	p = "HABCDS"
	6.times do |i|
		str.gsub!(i.to_s, p[i])
	end
	return str
end

#はじめの6つは mod 6、残りは mod 32 について記録している
@ransu = Array.new(100) { Array.new(6+32) { MiniSat::Var.new(@solver) } }

# 各ポケモンについてのありえる可能性（遺伝決定順序や破棄数）についての変数
@kanou = Array.new(poke.size) { Array.new }
kanou_comm = Array.new(poke.size) { Array.new }



# HABCDS のうちどれ？
def add_pos(pokenum, kouhonum, offset, value)
	if  !(0<= value && value < 6 ) then
		raise
	end
	
	@solver << [-@kanou[pokenum][kouhonum], @ransu[offset][value]]
	6.times do |v|
		if ( v != value ) then
			@solver << [-@kanou[pokenum][kouhonum], -@ransu[offset][v]]
		end
	end
	
	# mod 32のほうは偶奇の情報のみ
	add_kotaiguuki(pokenum, kouhonum, offset, value %2)
	
end


# HABCDS のうちどれ？（二通り考えられる場合）
def add_pos2(pokenum, kouhonum, offset, value1, value2)
	if  !(0<= value1 && value1 < 6 ) then
		raise
	end
	if  !(0<= value2 && value2 < 6 ) then
		raise
	end
	if ( value1 == value2 ) then
		raise
	end

	@solver << [-@kanou[pokenum][kouhonum], @ransu[offset][value1], @ransu[offset][value2]]
	6.times do |v|
		if ( v != value1 && v != value2 ) then
			@solver << [-@kanou[pokenum][kouhonum], -@ransu[offset][v]]
		end
	end
	
	
	# mod 32のほうは偶奇の情報のみ。二通りの偶奇が違うなら、情報がないことになるのでとばす
	if ( value1 %2 == value2 %2 ) then
		add_kotaiguuki(pokenum, kouhonum, offset, value1 %2)
	end
end


# どっち親からの遺伝？
def add_oya(pokenum, kouhonum, offset, value)

	if value.is_a? String then
		if value == '<' then
			value = 0
		elsif value == '>' then
			value = 1
		else
			raise
		end
	else
		raise
	end
	
	# mod 6 のほう
	add_mod6guuki(pokenum, kouhonum, offset, value)
	
	# mod 32 のほう
	add_kotaiguuki(pokenum, kouhonum, offset, value)
end
	
def add_mod6guuki(pokenum, kouhonum, offset, value)
	if ![0,1].include?(value) then
		raise
	end
	
	yep = [-@kanou[pokenum][kouhonum]]
	
	6.times do |v|
		if ( v%2 == value ) then
			yep.push ( @ransu[offset][v] )
		else
			@solver << [-@kanou[pokenum][kouhonum], -@ransu[offset][v]]
		end
	end
	
	@solver << yep
end


def add_kotai(pokenum, kouhonum, offset, value)

	if  !(0<= value && value < 32 ) then
		raise
	end
	
	32.times do |v|
		if ( v == value ) then
			@solver << [-@kanou[pokenum][kouhonum], @ransu[offset][6+v]]
		else
			@solver << [-@kanou[pokenum][kouhonum], -@ransu[offset][6+v]]
		end
	end
	
	# mod 6 のほう
	add_mod6guuki(pokenum, kouhonum, offset, value%2)
end
# mod 32 のほうの条件について、偶奇だけ入れる
def add_kotaiguuki(pokenum, kouhonum, offset, value)
	if [0,1].include?(value) then
	else
		raise
	end
	
	yep = [-@kanou[pokenum][kouhonum]]
	
	32.times do |v|
		if ( v%2 == value ) then
			yep.push ( @ransu[offset][6+v] )
		else
			@solver << [-@kanou[pokenum][kouhonum], -@ransu[offset][6+v]]
		end
	end
	
	@solver << yep
end


for pokenum in 0...(poke.size)

	poke_offset = pokenum #num回拒否したときの乱数スタート位置

	#遺伝した箇所の候補
	idenoya = ""
	6.times { |i|
		if ( poke[pokenum][i] == oya[0][i] ) then idenoya += "<"
		elsif ( poke[pokenum][i] == oya[1][i] ) then idenoya += ">"
		else
			idenoya += "."
		end
	}
	
	idens = (0...6).select { |i| ['<','>'].include?(idenoya[i]) }
	if idens.size < 3 then
		raise
	end
	
	# 遺伝の順番、（個体値一致が4つ以上ある場合の）遺伝箇所を列挙
	idens.permutation(3).each{ |iden|
	
		(0..MAX_SUTE1).each { |first_sute| #一個目の遺伝箇所を決定した後に破棄した回数
			(0..MAX_SUTE2).each { |second_sute| #二個目の遺伝箇所を決定した後に破棄した回数
				
				@kanou[pokenum].push(MiniSat::Var.new(@solver))
				kouhonum = @kanou[pokenum].size - 1
				
				
				teststring = ""
				
				offset = 0
				
				# 遺伝一個目
				add_pos(pokenum, kouhonum, poke_offset+0, iden[0])
				add_oya(pokenum, kouhonum, poke_offset+1, idenoya[iden[0]])
				offset += 2
				teststring += iden[0].to_s + idenoya[iden[0]].to_s
				
				# 一個目の後に捨てられる（＝一個目の遺伝箇所と同じ）
				first_sute.times do
					add_pos(pokenum, kouhonum, poke_offset+offset, iden[0])
					offset += 1
					teststring += "_"
				end
				
				# 遺伝二個目
				add_pos(pokenum, kouhonum, poke_offset+offset+0, iden[1])
				add_oya(pokenum, kouhonum, poke_offset+offset+1, idenoya[iden[1]])
				offset += 2
				teststring += iden[1].to_s + idenoya[iden[1]].to_s
				
				# 二個目の後に捨てられる
				second_sute.times do
					add_pos2(pokenum, kouhonum, poke_offset+offset, iden[0], iden[1])
					teststring += "|"
					offset += 1
				end
				
				# 遺伝三個目
				add_pos(pokenum, kouhonum, poke_offset+offset+0, iden[2])
				add_oya(pokenum, kouhonum, poke_offset+offset+1, idenoya[iden[2]])
				offset += 2
				teststring += iden[2].to_s + idenoya[iden[2]].to_s
				
				# 個体値の実数。遺伝済み箇所は飛ばす
				(0...6).each { |i|
					if !iden.include?(i) then
						add_kotai(pokenum, kouhonum, poke_offset+offset+i, poke[pokenum][i])
						teststring += "-+"[poke[pokenum][i]%2]
					else
						teststring += "."
					end
				}
				offset += 6
				if offset != teststring.size then
					raise
				end
				
				kanou_comm[pokenum].push({:text => pokenum.to_s + ": #{first_sute} #{second_sute} " + (" "*pokenum) + pretty(teststring),
				                          :sute => [first_sute, second_sute] } )
			}
		}
	
	}
	
	#ある孵化についての可能性のうち、ちょうどひとつだけ
	## https://www.cs.cmu.edu/~wklieber/papers/2007_efficient-cnf-encoding-for-selecting-1.pdf
	@solver << @kanou[pokenum]
	(@kanou[pokenum].size).times { |i|
		i.times { |j|
			@solver << [-@kanou[pokenum][i], -@kanou[pokenum][j]]
		}
	}
end


sute_kouho = -1
bit_kouho = Array.new(30, 3)

while true do

m = @solver.solve

p "--------"

if m then
	
	further = []
	
	poke.size.times { |num|
		@kanou[num].size.times { |k|
			if ( m[ @kanou[num][k] ] ) then
				p kanou_comm[num][k][:text]
				
				# 残りの可能性を探るために 条件を追加
				further.push( -@kanou[num][k] )
				
				
				if num == 0 then
					sutesum = kanou_comm[num][k][:sute].inject(:+)
					if sute_kouho == -1 then
						sute_kouho = sutesum
					elsif sute_kouho != sutesum then
						sute_kouho = -2
					end
				end
			end
		}
	}
	
	@solver << further
	
	str = "bitstream: "
	30.times { |i|
		may = 0
		if m[@ransu[i][0]] || m[@ransu[i][2]] || m[@ransu[i][4]] then
			may += 1
		end
		if m[@ransu[i][1]] || m[@ransu[i][3]] || m[@ransu[i][5]] then
			may += 2
		end
		
		if may == 3 then
			raise
		end
		
		str += "?01!"[may]
		
		bit_kouho[i] = bit_kouho[i] & may
	}
	
	p str
	
	
else

	if sute_kouho == -1 then
		p "Unsatisfied. Wrong IV?"
		exit
	end
	
	p "Finished"
	
	eff = 0
	str = "bit: "
	30.times { |i|
		str += "?01!"[bit_kouho[i]]
		if [1,2].include?(bit_kouho[i]) then
			eff += 1
		end
	}
	
	p str
	p "#{eff} bits"
	
	p "Haki Ransu in first Fuka: " + ((sute_kouho < 0) ? "Undetermined. Do more Rodeo!!" : sute_kouho.to_s)
	exit
end

end
